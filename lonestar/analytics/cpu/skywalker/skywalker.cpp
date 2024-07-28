
#include "skywalker.h"
#include "llvm/Support/CommandLine.h"

#include "sssp.hpp"

#include <iostream>
#include <vector>

static const char* name = "Skywalker";
static const char* desc = "";
static const char* url  = "skywalker";

void gen_partitions(idx_t* num_nodes, std::vector<idx_t>& parts,
                    std::vector<idx_t>& xadj, std::vector<idx_t>& adjncy) {

  idx_t options[METIS_NOPTIONS];

  idx_t ncon   = 1;
  idx_t nparts = NPARTS;
  idx_t objval;

  METIS_SetDefaultOptions(options);

  int res = METIS_PartGraphKway(num_nodes, &ncon, xadj.data(), adjncy.data(),
                                NULL, NULL, NULL, &nparts, NULL, NULL, options,
                                &objval, parts.data());

  if (res == METIS_OK) {
    galois::gInfo("METIS Partitioning Complete\n");
  } else {
    galois::gInfo("METIS Partitioning Failed\n");
  }
}

uint64_t sumResetLargeArray(VertexList<uint8_t>& arr) {
  galois::GAccumulator<uint64_t> sum;
  galois::do_all(
      galois::iterate(size_t{0}, arr.size()),
      [&](size_t i) {
        sum += arr.getData(i, galois::MethodFlag::READ);
        arr[i] = 0;
      },
      galois::steal(), galois::no_stats());

  return sum.reduce();
}

template <typename T>
GraphAlgorithm<T>::GraphAlgorithm(std::string& input) {
  galois::graphs::readGraph(graph, input);
  galois::gInfo("Read ", graph.size(), " nodes ", graph.sizeEdges(), " edges");

#ifdef METIS_SCHEME

  std::vector<idx_t> partitions(graph.size());
  std::vector<idx_t> xadj;
  std::vector<idx_t> adjncy;

  xadj.reserve(graph.size() + 1);

  idx_t curr_offset = 0;
  xadj.push_back(curr_offset);

  for (auto n : graph) {
    LNode<T>& data = graph.getData(n, galois::MethodFlag::UNPROTECTED);
    curr_offset += data.num_out_edges;
    xadj.push_back(curr_offset);

    for (auto e : graph.edges(n)) {
      adjncy.push_back(graph.getEdgeDst(e));
    }
  }

  idx_t num_nodes = graph.size();
  gen_partitions(&num_nodes, partitions, xadj, adjncy);
  adjncy.clear();
  xadj.clear();

#endif

  for (unsigned int i = 0; i < NPARTS; ++i) {
    ptn_sizes.construct(i, galois::GAccumulator<uint64_t>());
    ptn_update_vtxs.construct(i, VertexList<uint8_t>());
    ptn_update_vtxs[i].allocate(graph.size());
  }

  partition_ids.allocate(graph.size());
  curr_values.allocate(graph.size());
  prev_values.allocate(graph.size());
  agg_values.allocate(graph.size());
  out_degrees.allocate(graph.size());

  galois::do_all(
      galois::iterate(graph),
      [&](GNode<T> n) {
        curr_values[n] = 0;
        prev_values[n] = 0;
        agg_values[n]  = 0;
        out_degrees[n] =
            std::distance(graph.edges(n).begin(), graph.edges(n).end());

#ifdef METIS_SCHEME

        partition_ids[n] = partitions[n];
#else
        partition_ids[n] = n % NPARTS;
#endif
        ptn_sizes[partition_ids.getData(n, galois::MethodFlag::READ)] += 1;
      },
      galois::steal(), galois::loopname("Init Graph Strtucture"));

  for (unsigned int i = 0; i < NPARTS; ++i) {
    ptn_mirrors.construct(i, VertexList<T>());
    ptn_mirrors[i].allocate(graph.size());

    galois::do_all(
        galois::iterate(graph), [&](GNode<T> n) { ptn_mirrors[i][n] = 0; },
        galois::steal(), galois::loopname("Init Mirrors"));
  }

  dm_ndp_enabled.reset();
  dm_ndp_disabled.reset();

#ifdef METIS_SCHEME

  partitions.clear();

#endif
}

template <typename T>
void GraphAlgorithm<T>::run() {

  do {
    rounds++;

    // Update Mirrors
    galois::do_all(
        galois::iterate(frontier),
        [&](GNode<T>& src) {
          ptn_mirrors[partition_ids[src]][src] =
              curr_values.getData(src, galois::MethodFlag::READ);
          frontier_size += 1;

#ifdef SKYWALKER_DEBUG
          galois::gInfo("Updating mirror ", src, " with ",
                        curr_values.getData(src, galois::MethodFlag::READ));
#endif
        },
        galois::steal(), galois::loopname("Update Mirrors"));

#ifdef ITER_STATS
    galois::gInfo("Frontier Size: ", frontier_size.reduce());
#endif

    dm_ndp_disabled += frontier_size.reduce() * vtx_size;
    dm_ndp_enabled += frontier_size.reduce() * kv_size;

    ndp_disabled_iter += frontier_size.reduce() * vtx_size;
    ndp_enabled_iter += frontier_size.reduce() * kv_size;

    // Generate Updates
    for (unsigned int i = 0; i < NPARTS; ++i) {
      galois::do_all(
          galois::iterate(frontier),
          [&](GNode<T>& mirror) {
            if (partition_ids.getData(mirror, galois::MethodFlag::READ) != i) {
              return;
            }

            dm_ndp_disabled += std::distance(graph.edges(mirror).begin(),
                                             graph.edges(mirror).end()) *
                               vtx_size;

            ndp_disabled_iter += std::distance(graph.edges(mirror).begin(),
                                               graph.edges(mirror).end()) *
                                 vtx_size;

            generateUpdates(i, mirror);
          },
          galois::loopname("Generate Updates"));
    }

    frontier.clear();
    frontier_size.reset();

    // Aggregate Updates
    for (unsigned int i = 0; i < NPARTS; ++i) {
      galois::do_all(
          galois::iterate(size_t{0}, graph.size()),
          [&](size_t src) {
            if (!ptn_mirrors[i].isUpdated(src)) {
              return;
            }

            ptn_mirrors[i].resetUpdate(src);

            aggregateUpdates(
                agg_values[src],
                ptn_mirrors[i].getData(src, galois::MethodFlag::READ));

            ptn_update_vtxs[i][src] = 1;
          },
          galois::steal(), galois::loopname("Aggregate Updates"));

      uint64_t unique_vtxs = sumResetLargeArray(ptn_update_vtxs[i]);

#ifdef ITER_STATS
      galois::gInfo("Partition[", i, "] Unique Vtxs: ", unique_vtxs);
#endif

      dm_ndp_enabled += unique_vtxs * kv_size;
      ndp_enabled_iter += unique_vtxs * kv_size;
    }

    // Apply Updates
    galois::do_all(
        galois::iterate(graph),
        [&](GNode<T> src) {
          if (agg_values.getData(src, galois::MethodFlag::READ) == INT64_MAX) {
            return;
          }

          if (agg_values.getData(src, galois::MethodFlag::READ) <
              curr_values.getData(src, galois::MethodFlag::READ)) {
            curr_values[src] =
                agg_values.getData(src, galois::MethodFlag::READ);
            frontier.push(src);

#ifdef SKYWALKER_DEBUG
            galois::gInfo("Applying update to ", src, " with ",
                          agg_values.getData(src, galois::MethodFlag::READ));
#endif
          }

          agg_values[src] = INT64_MAX;
        },
        galois::steal(), galois::loopname("Apply Updates"));

#ifdef ITER_STATS
    galois::gInfo("Iter[", rounds,
                  "] NDP Enabled: ", ndp_enabled_iter.reduce());
    galois::gInfo("Iter[", rounds,
                  "] NDP Disabled: ", ndp_disabled_iter.reduce());
#endif

    ndp_enabled_iter.reset();
    ndp_disabled_iter.reset();

  } while (terminate() && rounds < MAX_ROUNDS);
}

template <typename T>
void GraphAlgorithm<T>::reportStats() {
  galois::reportPageAlloc("MeminfoPost");
  galois::runtime::reportStat_Single("Skywalker", "Iterations", rounds);
  galois::runtime::reportStat_Single(
      "Skywalker", "Data Movement (NDP Enabled) ", dm_ndp_enabled.reduce());
  galois::runtime::reportStat_Single(
      "Skywalker", "Data Movement (NDP Disabled)", dm_ndp_disabled.reduce());

#ifdef SKYWALKER_DEBUG
  // Print out the results
  for (auto n : graph) {
    std::cout << n << " " << curr_values[n] << "\n";
  }
#endif
}

int main(int argc, char** argv) {
  galois::SharedMemSys G;
  LonestarStart(argc, argv, name, desc, url, &inputFile);

  galois::StatTimer totalTime("TimerTotal");
  totalTime.start();

  galois::StatTimer execTime("Timer_0");
  execTime.start();

  SSSP<int64_t> algo_impl(inputFile);

  algo_impl.init();
  algo_impl.reportStats();

  algo_impl.run();
  algo_impl.reportStats();

  execTime.stop();
  totalTime.stop();

  return 0;
}
