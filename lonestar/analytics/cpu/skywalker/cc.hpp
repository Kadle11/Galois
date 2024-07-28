#ifndef SKYWALKER_CC_HPP
#define SKYWALKER_CC_HPP

#include "skywalker.h"

template <typename T>
class CC : public GraphAlgorithm<T> {
public:
  CC(std::string& input);

  void init();
  void generateUpdates(unsigned int& ptn_id, GNode<T>& mirror);
  void applyUpdates();
  void aggregateUpdates(T& value, const T& temp_val);
  void updateFrontier();
  bool terminate();
};

template <typename T>
CC<T>::CC(std::string& input) : GraphAlgorithm<T>(input) {}

template <typename T>
void CC<T>::init() {
  galois::do_all(
      galois::iterate(this->graph),
      [&](GNode<T> n) {
        this->curr_values[n] = n;
        this->prev_values[n] = n;
        this->agg_values[n]  = INT64_MAX;

        for (auto ii = 0; ii < NPARTS; ++ii) {
          this->ptn_updates[ii][n] = INT64_MAX;
        }

        this->frontier.push(n);
      },
      galois::steal(), galois::loopname("Init CC"));

}

template <typename T>
void CC<T>::generateUpdates(unsigned int& ptn_id, GNode<T>& mirror) {
  T new_dist = this->curr_values.getData(mirror);
  for (auto ii : this->graph.edges(mirror)) {
    GNode<T> dst = this->graph.getEdgeDst(ii);
    this->ptn_updates[ptn_id].minUpdate(dst, new_dist);
  }
}

template <typename T>
void CC<T>::applyUpdates() {
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
void CC<T>::aggregateUpdates(T& value, const T& temp_val) {
  value = std::min(value, temp_val);
}

template <typename T>
void CC<T>::updateFrontier() {}

template <typename T>
bool CC<T>::terminate() {
  return !this->frontier.empty();
}

// Explicit template instantiation
template class CC<uint64_t>;

#endif // SKYWALKER_CC_HPP
