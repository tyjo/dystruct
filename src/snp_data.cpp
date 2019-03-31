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

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <cassert>
#include <map>
#include <set>
#include <vector>
#include "snp_data.h"

#include <iostream>

using std::pair;
using std::map;
using std::vector;

SNPData::SNPData(const std_vector3<short> *snps, vector<int> sample_gen, double hold_out_proportion, int hold_out_seed)
    : sample_gen(sample_gen)
{
    // fixing the random seed will fix the hold out set across runs
    // so that the held-out log likelihood can be compared
    boost::random::mt19937 gen(hold_out_seed);
    this->snps = snps;

    // select hold_out_proportion SNP locations to put into a hold out set.
    // at each selected location, a genotype of a single individual is
    // held out.
    int total_individuals = 0;
    for (size_t t = 0; t < (*snps).size(); ++t)
        total_individuals += (*snps)[t].size();

    size_t nloci = (*snps)[0][0].size();
    int nlocations = nloci * hold_out_proportion;
    std::set<int> picked;
    boost::random::uniform_int_distribution<int> ldist(0, nloci - 1);
    boost::random::uniform_int_distribution<int> tdist(0, (*snps).size() - 1);

    int count = 0;
    for (int i = 0; i < nlocations; ++i) {

        // pick a SNP location
        int draw = ldist(gen);

        // make sure we don't pick the same location twice
        while (picked.find(draw) != picked.end()) {
            draw = ldist(gen);
        }
        picked.insert(draw);

        // pick a time
        int t = tdist(gen);;

        // pick an individual
        boost::random::uniform_int_distribution<int> idist(0, (*snps)[t].size() - 1);
        int individual = idist(gen);

        // make sure we have data
        while (missing(t, individual, draw)) {
            t = tdist(gen);
            boost::random::uniform_int_distribution<int> idist(0, (*snps)[t].size() - 1);
            individual = idist(gen);
        }

        if (ho.find(t) == ho.end()) {
            ho[t] = vector<vector<int> >();
            for (size_t d = 0; d < (*snps)[t].size(); ++d) {
                ho[t].push_back(vector<int>());
            }
        }
        ho[t][individual].push_back(draw);
        count += 1;
    }

    // check
    assert(count == nlocations);

    // check if individual is hemizygous
    for (size_t t = 0; t < (*snps).size(); ++t) {
        hemi.push_back(vector<bool>());
        for (size_t d = 0; d < (*snps)[t].size(); ++d) {
            bool is_hemizygous = true;
            for (size_t l = 0; l <  (*snps)[t][d].size(); ++l) {
                if ((*snps)[t][d][l] == 1) {
                    is_hemizygous = false;
                }
            }
            hemi[t].push_back(is_hemizygous);
        }
    }
}



int SNPData::max_individuals() const
{
    int max = 0;
    for (size_t t = 0; t < total_time_steps(); ++t) {
        if ((int)(*snps)[t].size() > max)
            max = (int)(*snps)[t].size();
    }
    return max;
}