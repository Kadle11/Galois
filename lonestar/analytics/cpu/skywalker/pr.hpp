#ifndef SKYWALKER_PR_HPP
#define SKYWALKER_PR_HPP

#include "skywalker.h"

#define TOLERANCE      1e-3 
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
};

template <typename T>
PR<T>::PR(std::string& input) : GraphAlgorithm<T>(input) {}

template <typename T>
void PR<T>::init() {
  galois::do_all(
      galois::iterate(this->graph),
      [&](GNode<T> n) {
        LNode<T>& data =
            this->graph.getData(n, galois::MethodFlag::UNPROTECTED);
        data.curr_val = 1.0 - DAMPING_FACTOR;
        data.prev_val = 1.0 - DAMPING_FACTOR;

        this->frontier.push(n);
      },
      galois::steal(), galois::loopname("Init PR"));
}

template <typename T>
void PR<T>::generateUpdates(unsigned int& ptn_id, GNode<T>& mirror) {
  LNode<T>& mdata = this->graph.getData(mirror, galois::MethodFlag::UNPROTECTED);
  
  for (auto ii : this->graph.edges(mirror)) {
    GNode<T> dst = this->graph.getEdgeDst(ii);

    this->ptn_updates[ptn_id].push(VertexUpdates(dst, update_val));
  }
}

#endif // SKYWALKER_PR_HPP