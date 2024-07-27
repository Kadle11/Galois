#ifndef SKYWALKER_PR_HPP
#define SKYWALKER_PR_HPP

#include "skywalker.h"

#define TOLERANCE 1e-3
#define DAMPING_FACTOR 0.85

template <typename T>
class PR : public GraphAlgorithm<T> {
public:
  PR(std::string& input);

  void init();
  void generateUpdates(unsigned int& ptn_id, GNode<T>& mirror);
  void applyUpdates();
  void aggregateUpdates(T& value, T& temp_val);
  void updateFrontier();
  bool terminate();

private:
  galois::LargeArray<T> pr_val;
  galois::LargeArray<T> deltas;
};

template <typename T>
PR<T>::PR(std::string& input) : GraphAlgorithm<T>(input) {}

template <typename T>
void PR<T>::init() {

  pr_val.allocateInterleaved(this->graph.size());
  deltas.allocateInterleaved(this->graph.size());

  galois::do_all(
      galois::iterate(this->graph),
      [&](GNode<T> n) {
        LNode<T>& data =
            this->graph.getData(n, galois::MethodFlag::UNPROTECTED);
        data.curr_val = (data.num_out_edges > 0) ? 1.0 - DAMPING_FACTOR : 0;
        deltas[n]     = DAMPING_FACTOR * data.curr_val / data.num_out_edges;

        pr_val[n] += data.curr_val;

        this->frontier.push(n);
      },
      galois::steal(), galois::loopname("Init PR"));
}

template <typename T>
void PR<T>::generateUpdates(unsigned int& ptn_id, GNode<T>& mirror) {

  for (auto ii : this->graph.edges(mirror)) {
    GNode<T> dst = this->graph.getEdgeDst(ii);
    this->ptn_updates[ptn_id].push(VertexUpdates(dst, deltas[mirror]));
  }
}

template <typename T>
void PR<T>::aggregateUpdates(T& value, T& temp_val) {
  value += temp_val;
}

template <typename T>
void PR<T>::applyUpdates() {
#ifdef SKYWALKER_DEBUG
  for (auto n : this->graph) {
    LNode<T>& data = this->graph.getData(n, galois::MethodFlag::UNPROTECTED);
    std::cout << n << " " << pr_val[n] << " " << data.curr_val << std::endl;
  }
#endif

  galois::do_all(
      galois::iterate(this->frontier),
      [&](GNode<T> n) {
        LNode<T>& data =
            this->graph.getData(n, galois::MethodFlag::UNPROTECTED);

        if (data.curr_val > TOLERANCE) {
          pr_val[n] += data.curr_val;
          deltas[n]     = DAMPING_FACTOR * data.curr_val / data.num_out_edges;
          data.prev_val = 0;
        }
      },
      galois::steal(), galois::loopname("Apply PR"));
}

template <typename T>
void PR<T>::updateFrontier() {
  galois::do_all(
      galois::iterate(this->graph),
      [&](GNode<T> n) {
        LNode<T>& data =
            this->graph.getData(n, galois::MethodFlag::UNPROTECTED);

        if (data.agg_val == 0) {
          return;
        }

        data.curr_val    = data.agg_val;
        T& prev_resiudal = data.prev_val;

        if (data.curr_val > TOLERANCE && prev_resiudal < TOLERANCE) {
          this->frontier.push(n);
          prev_resiudal = data.curr_val;
          data.agg_val  = 0;
        }
      },
      galois::steal(), galois::loopname("Update Frontier"));
}

template <typename T>
bool PR<T>::terminate() {

#ifdef SKYWALKER_DEBUG
  if (this->frontier.empty()) {
    for (auto n : this->graph) {
      std::cout << n << " " << pr_val[n] << std::endl;
    }
  }
#endif

  return !this->frontier.empty();
}

// Explicit template instantiation
template class PR<float>;
template class PR<double>;

#endif // SKYWALKER_PR_HPP