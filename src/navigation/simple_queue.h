// Copyright (c) 2018 joydeepb@cs.umass.edu
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <algorithm>
#include <cstddef>
#include <deque>
#include <sys/types.h>
#include <utility>
#include <vector>
#include <tuple>

#ifndef MUTABLE_QUEUE
#define MUTABLE_QUEUE

using std::deque;
using std::pair;
using std::make_pair;
using std::tuple;
using std::make_tuple;

template<class CycleNum, class Velocity, class Curvature>
class SimpleQueue {
 private:
  public:
  // Insert a new value, with the specified priority. If the value
  // already exists, its priority is updated.
  void Push(const CycleNum& n, const Velocity& v, const Curvature& c) {
    values_.push_back(make_tuple(n, v, c));
  }
  
  // Retreive the value with the highest priority.
  pair<Velocity, Curvature> Pop() {
    if (values_.empty()) {
      fprintf(stderr, "ERROR: Pop() called on an empty queue!\n");
      exit(1);
    }
    const Velocity v = std::get<1>(values_.front());
    const Curvature p = std::get<2>(values_.front());
    values_.pop_front();
    return make_pair(v, p);
  }

  // Returns true iff the priority queue is empty.
  bool Empty() {
    return values_.empty();
  }

  // Returns true iff the provided value is already on the queue.
  bool Exists(const CycleNum& n) {
    for (const auto& x : values_) {
      if (x.first == n) return true;
    }
    return false;
  }

  unsigned Size() {
    return (unsigned) values_.size();
  }

  deque<tuple<CycleNum, Velocity, Curvature>> values_;
};

#endif  // MUTABLE_QUEUE
