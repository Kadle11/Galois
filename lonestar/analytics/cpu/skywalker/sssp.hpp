#ifndef SKYWALKER_SSSP_HPP
#define SKYWALKER_SSSP_HPP

#include "skywalker.h"

template <typename T>
class SSSP : public GraphAlgorithm<T> {
public:
  SSSP(std::string& input);

  void init();
  void generateUpdates(unsigned int& ptn_id, GNode<T>& mirror);
  void applyUpdates();
  void aggregateUpdates(T& value, const T& temp_val);
  void updateFrontier();
  bool terminate();
};

template <typename T>
SSSP<T>::SSSP(std::string& input) : GraphAlgorithm<T>(input) {}

template <typename T>
void SSSP<T>::init() {
  galois::do_all(
      galois::iterate(this->graph),
      [&](GNode<T> n) {
        this->curr_values[n] = INT64_MAX;
        this->prev_values[n] = INT64_MAX;
        this->agg_values[n]  = INT64_MAX;

        for (auto ii = 0; ii < NPARTS; ++ii) {
          this->ptn_updates[ii][n] = INT64_MAX;
        }
      },
      galois::steal(), galois::loopname("Init SSSP"));

  this->frontier.push(0);
  this->curr_values[0] = 0;
}

template <typename T>
void SSSP<T>::generateUpdates(unsigned int& ptn_id, GNode<T>& mirror) {
  T new_dist = this->curr_values.getData(mirror) + 1;
  for (auto ii : this->graph.edges(mirror)) {
    GNode<T> dst = this->graph.getEdgeDst(ii);
    this->ptn_updates[ptn_id].minUpdate(dst, new_dist);
  }
}

template <typename T>
void SSSP<T>::applyUpdates() {
  // Apply Updates
  galois::do_all(
      galois::iterate(this->graph),
      [&](GNode<T> src) {
        if (this->agg_values.getData(src) == INT64_MAX) {
          return;
        }

        if (this->agg_values.getData(src) < this->curr_values.getData(src)) {
          this->curr_values[src] = this->agg_values.getData(src);
          this->frontier.push(src);

#ifdef SKYWALKER_DEBUG
          galois::gInfo("Applying update to ", src, " with ",
                        this->agg_values.getData(src));
#endif
        }

        this->agg_values[src] = INT64_MAX;
      },
      galois::steal(), galois::loopname("Apply Updates"));
}

template <typename T>
void SSSP<T>::aggregateUpdates(T& value, const T& temp_val) {
  value = std::min(value, temp_val);
}

template <typename T>
void SSSP<T>::updateFrontier() {}

template <typename T>
bool SSSP<T>::terminate() {
  return !this->frontier.empty();
}

// Explicit template instantiation
template class SSSP<uint64_t>;

#endif // SKYWALKER_SSSP_HPP
