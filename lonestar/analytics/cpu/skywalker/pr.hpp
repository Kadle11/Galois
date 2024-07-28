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
  void aggregateUpdates(T& value, const T& temp_val);
  void updateFrontier();
  bool terminate();

private:
  VertexList<T> pr_val;
  VertexList<T> deltas;
};

template <typename T>
PR<T>::PR(std::string& input) : GraphAlgorithm<T>(input) {}

template <typename T>
void PR<T>::init() {

  pr_val.allocate(this->graph.size());
  deltas.allocate(this->graph.size());

  galois::do_all(
      galois::iterate(this->graph),
      [&](GNode<T> n) {
        this->curr_values[n] =
            (this->out_degrees.getData(n) > 0) ? 1.0 - DAMPING_FACTOR : 0;

        deltas[n] = DAMPING_FACTOR * this->curr_values.getData(n) /
                    this->out_degrees.getData(n);

        pr_val[n] = this->curr_values.getData(n);

        this->curr_values[n] = 0;
        this->frontier.push(n);
      },
      galois::steal(), galois::loopname("Init PR"));
}

template <typename T>
void PR<T>::generateUpdates(unsigned int& ptn_id, GNode<T>& mirror) {
  T delta = deltas.getData(mirror);

  for (auto ii : this->graph.edges(mirror)) {
    GNode<T> dst = this->graph.getEdgeDst(ii);

#ifdef SKYWALKER_DEBUG
    galois::gInfo("Updating mirror of ", dst, " with ",
                  this->ptn_updates[ptn_id].getData(dst), " + ", delta);
#endif

    this->ptn_updates[ptn_id].addUpdate(dst, delta);
  }
}

template <typename T>
void PR<T>::aggregateUpdates(T& value, const T& temp_val) {
#ifdef SKYWALKER_DEBUG
  galois::gInfo("Aggregating ", value, " with ", temp_val);
#endif
  value += temp_val;
}

template <typename T>
void PR<T>::applyUpdates() {
#ifdef SKYWALKER_DEBUG
  for (auto n : this->graph) {
    std::cout << n << " " << pr_val[n] << " " << this->curr_values.getData(n)
              << std::endl;
  }
#endif

  galois::do_all(
      galois::iterate(this->frontier),
      [&](GNode<T> n) {
        if (this->curr_values.getData(n) > TOLERANCE) {
          pr_val.addUpdate(n, this->curr_values.getData(n));
          deltas[n] = DAMPING_FACTOR * this->curr_values.getData(n) /
                      this->out_degrees.getData(n);

          this->prev_values[n] = 0;
          this->curr_values[n] = 0;
        }
      },
      galois::steal(), galois::loopname("Apply PR"));
}

template <typename T>
void PR<T>::updateFrontier() {
  galois::do_all(
      galois::iterate(this->graph),
      [&](GNode<T> n) {
        if (this->agg_values.getData(n) == 0) {
          return;
        }
#ifdef SKYWALKER_DEBUG

        galois::gInfo("Updating frontier vertex ", n, " with ",
                      this->agg_values.getData(n));
#endif
        this->curr_values.addUpdate(n, this->agg_values.getData(n));

        if (this->curr_values.getData(n) > TOLERANCE &&
            this->prev_values.getData(n) < TOLERANCE) {

          this->frontier.push(n);
          this->prev_values[n] = this->curr_values.getData(n);
        }

        this->agg_values[n] = 0;
        for (int i = 0; i < NPARTS; i++) {
          this->ptn_updates[i][n] = 0;
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

  return !this->frontier.empty() && this->rounds < 20;
}

// Explicit template instantiation
template class PR<float>;
template class PR<double>;

#endif // SKYWALKER_PR_HPP