#ifndef SKYWALKER_H
#define SKYWALKER_H

#include "galois/Galois.h"
#include "galois/AtomicHelpers.h"
#include "galois/Reduction.h"
#include "galois/PriorityQueue.h"
#include "galois/Timer.h"
#include "galois/graphs/LCGraph.h"
#include "galois/graphs/TypeTraits.h"
#include "Lonestar/BoilerPlate.h"

#include <iostream>
#include <vector>
#include <string>
#include <shared_mutex>

#include "metis.h"

// #define SKYWALKER_DEBUG
// #define ITER_STATS
// #define METIS_SCHEME

// #define NPARTS 2

namespace cll = llvm::cl;
static cll::opt<std::string>
    inputFile(cll::Positional, cll::desc("<input file>"), cll::Required);

template <typename T>
struct LNode {
  T curr_val;
  T prev_val;

  T agg_val;
  uint64_t num_out_edges;
};

template <typename T>
using Graph = typename galois::graphs::LC_CSR_Graph<
    uint64_t, uint32_t>::with_no_lockable<true>::type;

template <typename T>
using GNode = typename Graph<T>::GraphNode;

template <typename T>
struct dataElement : public galois::runtime::Lockable {

public:
  using reference = T&;

  T v;
  bool updated;

  reference getData() { return v; }
};

template <typename T>
struct VertexList {

public:
  VertexList() = default;

  typename dataElement<T>::reference operator[](const GNode<T>& n) {

    acquireNode(n, galois::MethodFlag::WRITE);
    return data[n].getData();
  }

  typename dataElement<T>::reference
  getData(const GNode<T>& n,
          galois::MethodFlag mflag = galois::MethodFlag::READ) {

    acquireNode(n, mflag);
    return data[n].getData();
  }

  void minUpdate(const GNode<T>& n, const T& val) {
    acquireNode(n);
    data[n].v       = std::min(data[n].v, val);
    data[n].updated = true;
  }

  void maxUpdate(const GNode<T>& n, const T& val) {
    acquireNode(n);
    data[n].v       = std::max(data[n].v, val);
    data[n].updated = true;
  }

  void addUpdate(const GNode<T>& n, const T& val) {
    acquireNode(n);
    data[n].v += val;
    data[n].updated = true;
  }

  void setUpdate(const GNode<T>& n, const T& val) {
    acquireNode(n);
    data[n].v       = val;
    data[n].updated = true;
  }

  void resetUpdate(const GNode<T>& n) {
    acquireNode(n);
    data[n].updated = false;
  }

  bool isUpdated(const GNode<T>& n) {
    acquireNode(n, galois::MethodFlag::READ);
    return data[n].updated;
  }

  void allocate(size_t size) { data.allocateInterleaved(size); }

  using iterator = boost::counting_iterator<GNode<T>>;
  iterator begin() { return iterator(0); }
  iterator end() { return iterator(data.size()); }

  size_t size() { return data.size(); }

private:
  galois::LargeArray<dataElement<T>> data;
  void acquireNode(const GNode<T>& node,
                   galois::MethodFlag mflag = galois::MethodFlag::WRITE) {

    assert(n < data.size());
    galois::runtime::acquire(&data[node], mflag);
  }
};

template <typename T>
class GraphAlgorithm {

public:
  GraphAlgorithm(std::string& input);

  void run();
  void reportStats();

  virtual void init()                                                  = 0;
  virtual void generateUpdates(unsigned int& ptn_id, GNode<T>& mirror) = 0;
  virtual void applyUpdates()                                          = 0;
  virtual void aggregateUpdates(T& value, const T& temp_val)           = 0;
  virtual void updateFrontier()                                        = 0;
  virtual bool terminate()                                             = 0;

protected:
  Graph<T> graph;

  galois::InsertBag<GNode<T>> frontier;
  galois::LazyArray<VertexList<T>, NPARTS> ptn_updates;
  VertexList<T> partition_ids;
  VertexList<T> agg_values;
  VertexList<T> curr_values;
  VertexList<T> prev_values;
  VertexList<T> out_degrees;

  // Structures for Telemetry
  galois::GAccumulator<uint64_t> dm_ndp_enabled;
  galois::GAccumulator<uint64_t> dm_ndp_disabled;
  galois::GAccumulator<uint64_t> frontier_size;
  galois::GAccumulator<uint64_t> ndp_enabled_iter;
  galois::GAccumulator<uint64_t> ndp_disabled_iter;
  galois::LazyArray<galois::GAccumulator<uint64_t>, NPARTS> ptn_sizes;
  galois::LazyArray<VertexList<uint8_t>, NPARTS> ptn_update_vtxs;

  uint64_t rounds                      = 0;
  static constexpr uint64_t vtx_size   = 8;
  static constexpr uint64_t kv_size    = 16;
  static const unsigned int MAX_ROUNDS = 100;
};

#endif // SKYWALKER_H
