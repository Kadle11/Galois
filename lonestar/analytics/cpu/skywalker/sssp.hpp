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
  void aggregateUpdates(T& value, T& temp_val);
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
        LNode<T>& data =
            this->graph.getData(n, galois::MethodFlag::UNPROTECTED);

        data.curr_val = INT64_MAX;
        data.prev_val = INT64_MAX;
        data.agg_val  = INT64_MAX;
      },
      galois::steal(), galois::loopname("Init SSSP"));

  this->frontier.push(0);
  this->graph.getData(0).curr_val = 0;
}

template <typename T>
void SSSP<T>::generateUpdates(unsigned int& ptn_id, GNode<T>& mirror) {
  for (auto ii : this->graph.edges(mirror)) {
    GNode<T> dst = this->graph.getEdgeDst(ii);

    T& curr_dist = this->ptn_mirrors[ptn_id][dst];
    T new_dist   = this->ptn_mirrors[ptn_id][mirror] + 1;

    if (new_dist < curr_dist) {
      this->ptn_updates[ptn_id].push(VertexUpdates(dst, new_dist));
      // this->ptn_mirrors[ptn_id][dst] = new_dist;
    }
  }
}

template <typename T>
void SSSP<T>::applyUpdates() {
  galois::do_all(
      galois::iterate(this->graph),
      [&](GNode<T> n) {
        LNode<T>& data =
            this->graph.getData(n, galois::MethodFlag::UNPROTECTED);

        if (data.curr_val > data.agg_val) {
          data.curr_val = data.agg_val;
        }
      },
      galois::steal(), galois::loopname("Apply Updates"));
}

template <typename T>
void SSSP<T>::aggregateUpdates(T& value, T& temp_val) {
  value = std::min(value, temp_val);
}

template <typename T>
void SSSP<T>::updateFrontier() {
  galois::do_all(
      galois::iterate(this->graph),
      [&](GNode<T> n) {
        LNode<T>& data =
            this->graph.getData(n, galois::MethodFlag::UNPROTECTED);
        if (data.curr_val < data.agg_val) {
          this->frontier.push(n);
        }
      },
      galois::steal(), galois::loopname("Update Frontier"));
}

template <typename T>
bool SSSP<T>::terminate() {
  return !this->frontier.empty();
}

// Explicit template instantiation
template class SSSP<uint64_t>;

#endif // SKYWALKER_SSSP_HPP
