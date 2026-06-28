// tree.h

#pragma once

#include <cstdint>

template <typename T, unsigned int D>
struct Tree {
  bool _isLeaf;
  // refinement _level, non-refinement _level = 0
  std::uint8_t _level;
  // index in block
  std::size_t _id;

  Tree<T, D>* _parent;
  Tree<T, D>** _children;

  Tree(std::size_t id, std::uint8_t level, Tree<T, D>* parent = nullptr, bool isLeaf = true)
      : _id(id), _level(level), _parent(parent), _isLeaf(isLeaf) {}
  ~Tree() {
    if (!_isLeaf) {
      for (unsigned int i = 0; i < D; ++i) {
        delete _children[i];
      }
      delete[] _children;
    }
  }
  // void refine() {
  //   if (_isLeaf) {
  //     _isLeaf = false;
  //     _children = new Tree<T, D>*[D];
  //     for (unsigned int i = 0; i < D; ++i) {
  //       _children[i] = new Tree<T, D>(D * _id + i, _level + 1, this);
  //     }
  //   }
  // }
};