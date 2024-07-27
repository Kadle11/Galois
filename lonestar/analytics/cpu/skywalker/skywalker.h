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

#define NPARTS 128

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
    LNode<T>, uint32_t>::with_no_lockable<true>::type;

template <typename T>
using GNode = typename Graph<T>::GraphNode;

template <typename T>
struct VertexMirror {
  T value;
  bool updated;

  std::shared_mutex _mtx;

  VertexMirror() : value(0), updated(false) {}
  T getValue() {
    std::shared_lock<std::shared_mutex> lock(_mtx);
    return value;
  }

  void initialize(const T& val) {
    std::lock_guard<std::shared_mutex> lock(_mtx);
    value = val;
  }

  void setValue(const T& val) {
    std::lock_guard<std::shared_mutex> lock(_mtx);
    value   = val;
    updated = true;
  }

  void minUpdate(const T& val) {
    std::lock_guard<std::shared_mutex> lock(_mtx);
    if (val < value) {
      value   = val;
      updated = true;
    }
  }

  void maxUpdate(const T& val) {
    std::lock_guard<std::shared_mutex> lock(_mtx);
    if (val > value) {
      value   = val;
      updated = true;
    }
  }

  void addUpdate(const T& val) {
    std::lock_guard<std::shared_mutex> lock(_mtx);
    value += val;
    updated = true;
  }

  bool isUpdated() {
    std::shared_lock<std::shared_mutex> lock(_mtx);
    return updated;
  }

  void reset() {
    std::lock_guard<std::shared_mutex> lock(_mtx);
    updated = false;
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
  galois::LazyArray<galois::LargeArray<VertexMirror<T>>, NPARTS> ptn_mirrors;
  galois::LargeArray<unsigned int> partition_ids;

  // Structures for Telemetry
  galois::GAccumulator<uint64_t> dm_ndp_enabled;
  galois::GAccumulator<uint64_t> dm_ndp_disabled;
  galois::GAccumulator<uint64_t> frontier_size;
  galois::GAccumulator<uint64_t> ndp_enabled_iter;
  galois::GAccumulator<uint64_t> ndp_disabled_iter;
  galois::LazyArray<galois::GAccumulator<uint64_t>, NPARTS> ptn_sizes;
  galois::LazyArray<galois::LargeArray<uint8_t>, NPARTS> ptn_update_vtxs;

  uint64_t rounds                      = 0;
  static constexpr uint64_t vtx_size   = 8;
  static constexpr uint64_t kv_size    = 16;
  static const unsigned int MAX_ROUNDS = 100;
};

#endif // SKYWALKER_H
