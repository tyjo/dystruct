/*
Copyright (C) 2017-2018 Tyler Joseph <tjoseph@cs.columbia.edu>

This file is part of Dystruct.

Dystruct is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Dystruct is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Dystruct.  If not, see <http://www.gnu.org/licenses/>.
*/

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