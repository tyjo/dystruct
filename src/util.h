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

#ifndef UTIL_DYSTRUCT_H
#define UTIL_DYSTRUCT_H

#include <string>
#include <vector>

#include "snp_data.h"
#include "vector_types.h"

void read_snp_matrix(std::string fname, std_vector3<short> *snps, std::vector<int>& gen_sampled, int nloci);
vector2<int> read_pop_labels(std::string fname, SNPData& snp_data);

#endif