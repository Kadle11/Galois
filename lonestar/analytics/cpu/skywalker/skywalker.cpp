#include "galois/Galois.h"
#include "galois/AtomicHelpers.h"
#include "galois/Reduction.h"
#include "galois/PriorityQueue.h"
#include "galois/Timer.h"
#include "galois/graphs/LCGraph.h"
#include "galois/graphs/TypeTraits.h"
#include "Lonestar/BoilerPlate.h"

#include "llvm/Support/CommandLine.h"

#include <iostream>
#include <vector>

// #define SKYWALKER_DEBUG

namespace cll = llvm::cl;

static const char* name              = "Skywalker";
static const char* desc              = "";
static const char* url               = "skywalker";
static const unsigned int NPARTS     = 4;
static const unsigned int MAX_ROUNDS = 100;

static cll::opt<std::string>
    inputFile(cll::Positional, cll::desc("<input file>"), cll::Required);

struct LNode {
  uint64_t curr_val;
  uint64_t prev_val;

  uint64_t agg_val;
};

using Graph =
    galois::graphs::LC_CSR_Graph<LNode, uint32_t>::with_no_lockable<true>::type;

using GNode = Graph::GraphNode;

struct VertexUpdates {
  GNode src;
  uint64_t val;
  VertexUpdates(const GNode& N, uint64_t W) : src(N), val(W) {}
  VertexUpdates() : src(), val(INT64_MAX) {}
};

int main(int argc, char** argv) {
  galois::SharedMemSys G;
  LonestarStart(argc, argv, name, desc, url, &inputFile);

  galois::StatTimer totalTime("TimerTotal");
  totalTime.start();

  Graph graph;

  std::cout << "Reading from file: " << inputFile << "\n";
  galois::graphs::readGraph(graph, inputFile);
  std::cout << "Nodes: " << graph.size() << ", Edges: " << graph.sizeEdges()
            << std::endl;

  galois::LargeArray<unsigned int> partition_ids;
  partition_ids.allocateInterleaved(graph.size());

  galois::InsertBag<GNode> frontier;
  galois::LazyArray<galois::InsertBag<VertexUpdates>, NPARTS> ptn_updates;
  galois::LazyArray<galois::LargeArray<uint64_t>, NPARTS> ptn_mirrors;

  // Structures for Telemetry
  galois::GAccumulator<uint64_t> dm_ndp_enabled;
  galois::GAccumulator<uint64_t> dm_ndp_disabled;
  galois::GAccumulator<uint64_t> frontier_size;
  galois::LazyArray<galois::GAccumulator<uint64_t>, NPARTS> ptn_sizes;

  uint64_t rounds             = 0;
  constexpr uint64_t vtx_size = 8;
  constexpr uint64_t kv_size  = 16;

  // Allocate LazyArrays
  for (unsigned int i = 0; i < NPARTS; ++i) {
    ptn_updates.construct(i, galois::InsertBag<VertexUpdates>());
    ptn_sizes.construct(i, galois::GAccumulator<uint64_t>());
  }

  galois::StatTimer execTime("Timer_0");
  execTime.start();

  galois::do_all(
      galois::iterate(graph),
      [&](GNode n) {
        LNode& data   = graph.getData(n, galois::MethodFlag::UNPROTECTED);
        data.curr_val = INT64_MAX;
        data.agg_val  = INT64_MAX;
        partition_ids.constructAt(n, n % NPARTS);
        ptn_sizes[partition_ids[n]] += 1;
      },
      galois::steal(), galois::loopname("Init Graph Strtucture"));

  for (unsigned int i = 0; i < NPARTS; ++i) {
    ptn_mirrors.construct(i, galois::LargeArray<uint64_t>());
    ptn_mirrors[i].allocateInterleaved(graph.size());
    galois::do_all(
        galois::iterate(graph),
        [&](GNode n) { ptn_mirrors[i].constructAt(n, INT64_MAX); },
        galois::steal(), galois::loopname("Init Mirrors"));
  }

  frontier.push(0);
  graph.getData(0).curr_val = 0;

  do {

    rounds++;

    // Update Mirrors
    galois::do_all(
        galois::iterate(frontier),
        [&](GNode& src) {
          LNode& sdata = graph.getData(src, galois::MethodFlag::UNPROTECTED);
          uint64_t val = sdata.curr_val;

          ptn_mirrors[partition_ids[src]][src] = val;
          frontier_size += 1;

#ifdef SKYWALKER_DEBUG
          galois::gInfo("Updating mirror ", src, " with ", val);
#endif
        },
        galois::steal(), galois::loopname("Update Mirrors"));

    dm_ndp_disabled += frontier_size.reduce() * vtx_size;
    dm_ndp_enabled += frontier_size.reduce() * kv_size;

    galois::gInfo("Frontier Size: ", frontier_size.reduce());

    // Generate Updates
    for (unsigned int i = 0; i < NPARTS; ++i) {
      galois::do_all(
          galois::iterate(frontier),
          [&](GNode& mirror) {
            if (partition_ids[mirror] != i) {
              return;
            }

            dm_ndp_disabled += std::distance(graph.edges(mirror).begin(),
                                             graph.edges(mirror).end()) *
                               vtx_size;

            for (auto ii : graph.edges(mirror)) {
              GNode dst = graph.getEdgeDst(ii);

              uint64_t& curr_dist = ptn_mirrors[i][dst];
              uint64_t new_dist   = ptn_mirrors[i][mirror] + 1;

              if (new_dist < curr_dist) {
                ptn_updates[i].push(VertexUpdates(dst, new_dist));
#ifdef SKYWALKER_DEBUG
                galois::gInfo("Partition[", i, "] Pushing update from ", mirror,
                              " to ", dst, " with ", new_dist);
#endif
              }
            }
          },
          galois::steal(), galois::loopname("Generate Updates"));
    }

    frontier.clear();
    frontier_size.reset();

    // Aggregate Updates
    for (unsigned int i = 0; i < NPARTS; ++i) {
      galois::do_all(
          galois::iterate(ptn_updates[i]),
          [&](VertexUpdates& ptn_update) {
            LNode& sdata =
                graph.getData(ptn_update.src, galois::MethodFlag::UNPROTECTED);
            sdata.agg_val = std::min(sdata.agg_val, ptn_update.val);
          },
          galois::steal(), galois::loopname("Aggregate Updates"));
      ptn_updates[i].clear();
    }

    // Apply Updates
    galois::do_all(
        galois::iterate(graph),
        [&](GNode src) {
          LNode& sdata = graph.getData(src, galois::MethodFlag::UNPROTECTED);
          if (sdata.agg_val == INT64_MAX) {
            return;
          }

          dm_ndp_enabled += kv_size;

          if (sdata.agg_val < sdata.curr_val) {
            sdata.curr_val = sdata.agg_val;
            frontier.push(src);

#ifdef SKYWALKER_DEBUG
            galois::gInfo("Applying update to ", src, " with ", sdata.agg_val);
#endif
          }

          sdata.agg_val = INT64_MAX;
        },
        galois::steal(), galois::loopname("Apply Updates"));

    // aggregate_updates.clear();

  } while (!frontier.empty() && rounds < MAX_ROUNDS);

  // Print out the results

#ifdef SKYWALKER_DEBUG
  for (auto n : graph) {
    std::cout << n << " " << graph.getData(n).curr_val << "\n";
  }
#endif

  execTime.stop();
  totalTime.stop();

  galois::reportPageAlloc("MeminfoPost");
  galois::runtime::reportStat_Single("Skywalker", "Iterations", rounds);
  galois::runtime::reportStat_Single(
      "Skywalker", "Data Movement (NDP Enabled) ", dm_ndp_enabled.reduce());
  galois::runtime::reportStat_Single(
      "Skywalker", "Data Movement (NDP Disabled)", dm_ndp_disabled.reduce());

  return 0;
}
