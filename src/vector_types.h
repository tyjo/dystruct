#ifndef VECTOR_TYPES_H
#define VECTOR_TYPES_H

#include <boost/multi_array.hpp>

template <typename T>
using vector2 = boost::multi_array<T, 2>;

template <typename T>
using vector3 = boost::multi_array<T, 3>;

template <typename T>
using vector4 = boost::multi_array<T, 4>;

template <typename T>
using std_vector3 = std::vector<vector2<T> >;

#endif