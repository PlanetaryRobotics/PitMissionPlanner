#pragma once
#ifndef PLANRANGER_PRIORITY_QUEUE_H
#define PLANRANGER_PRIORITY_QUEUE_H

#include <algorithm>
#include <cassert>
#include <functional>
#include <queue>
#include <utility>
#include <vector>

#include "flat_hash_map.hpp"

template <class VALUE, class PRIORITY>
class PriorityQueue {
  struct HeapEntry;

 public:
  VALUE top() { return vec[0].value; }
  PRIORITY top_priority() { return vec[0].priority; }

  void pop() { remove(top()); }

  void insert(const VALUE& v, const PRIORITY& p) {
    vec.push_back(HeapEntry{v, p});
    std::push_heap(vec.begin(), vec.end());
  }

  void remove(const VALUE& v) {
    // TRY NOT TO USE THIS. IT IS SLOW AS HELL.
    int i = 0;
    for (; i < vec.size(); ++i) {
      if (vec[i].value == v) {
        break;
      }
    }
    if (i == vec.size()) {
      return;
    }
    vec.erase(vec.begin() + i);
    std::make_heap(vec.begin(), vec.end());
  }

  bool empty() const { return vec.size() == 0; }
  std::size_t size() const { return vec.size(); }

 public:
  std::vector<HeapEntry> vec;

 private:
  struct HeapEntry {
    VALUE value;
    PRIORITY priority;

    bool operator<(const HeapEntry& rhs) const {
      // std::push_heap, std::pop_heap, etc.
      // construct a max. heap. I need a min. heap,
      // so this inequality is backwards from what you might expect.
      return priority > rhs.priority;
    }
  };
};

#endif  // PLANRANGER_PRIORITY_QUEUE_H
