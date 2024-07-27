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
  void aggregateUpdates(T& value, T& temp_val);
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
        LNode<T>& data =
            this->graph.getData(n, galois::MethodFlag::UNPROTECTED);

        data.curr_val = n;
        data.prev_val = n;
        data.agg_val  = INT64_MAX;

        this->frontier.push(n);

        // galois::gInfo("Node ", n, " has degree ", data.num_out_edges);
      },
      galois::steal(), galois::loopname("Init CC"));
}

template <typename T>
void CC<T>::generateUpdates(unsigned int& ptn_id, GNode<T>& mirror) {
  T& src_val = this->graph.getData(mirror).curr_val;
  for (auto ii : this->graph.edges(mirror)) {
    GNode<T> dst = this->graph.getEdgeDst(ii);

    T& curr_val = this->ptn_mirrors[ptn_id][dst];

    if (src_val < curr_val) {
      this->ptn_updates[ptn_id].push(VertexUpdates(dst, src_val));
      this->ptn_mirrors[ptn_id][dst] = src_val;
    }
  }
}

template <typename T>
void CC<T>::applyUpdates() {
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
void CC<T>::aggregateUpdates(T& value, T& temp_val) {
  value = std::min(value, temp_val);
}

template <typename T>
void CC<T>::updateFrontier() {
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
bool CC<T>::terminate() {
  return !this->frontier.empty();
}

// Explicit template instantiation
template class CC<uint64_t>;

#endif // SKYWALKER_CC_HPP
