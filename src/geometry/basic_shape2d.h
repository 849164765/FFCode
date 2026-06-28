#pragma once

#include "geometry/basic_geometry.h"

template <typename T>
class Circle final : public AABB<T, 2> {
 protected:
  T _Radius;

 public:
  Circle(T radius, const Vector<T, 2>& centre)
      : AABB<T, 2>(centre, 2 * radius, 2 * radius), _Radius(radius) {}
  bool isInside(const Vector<T, 2>& pt) const override {
    constexpr T eps = std::numeric_limits<T>::epsilon();
    return (pt - this->_center).getnorm2() <= _Radius * _Radius + eps;
  }
  template <typename S>
  bool IsInside(const Vector<S, 2>& pt) const {
    T eps = std::numeric_limits<decltype(T{} * S{})>::epsilon();
    return (pt - this->_center).getnorm2() <= _Radius * _Radius + eps;
  }
};