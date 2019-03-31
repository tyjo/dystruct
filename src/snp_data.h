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

#ifndef SNP_DATA_H
#define SNP_DATA_H

#include <algorithm>
#include <boost/random/mersenne_twister.hpp>
#include <map>
#include <utility>
#include <vector>

#include "vector_types.h"

#include <iostream>


class SNPData
{
    public:
        SNPData() { };
        SNPData(SNPData& d) : ho(d.ho), sample_gen(d.sample_gen), hemi(d.hemi) { this->snps = d.snps; }
        SNPData(const std_vector3<short> *snps, std::vector<int> sample_gen, double hold_out_proportion, int hold_out_seed);
        size_t total_time_steps() const                                   { return (*snps).size(); }
        size_t total_individuals(size_t time) const                       { return (*snps)[time].size(); }
        // total loci for individual indiv sampled at time t
        size_t total_loci(size_t time, size_t indiv) const                { return (*snps)[time][indiv].size(); }
        int    get_sample_gen(size_t time) const                          { return sample_gen[time]; }
        int    max_individuals() const;

        double genotype(size_t time, size_t indiv, size_t locus) const
        { 
            return (double)(*snps)[time][indiv][locus];
        }

        bool   missing(size_t time, size_t indiv, size_t locus) const
        {
            return (*snps)[time][indiv][locus] == (short)9;
        }

        bool   hold_out(int time, size_t indiv, int locus) const
        {
            if (ho.find(time) == ho.end())
                return false;
            else 
                return (find(ho.at(time)[indiv].begin(), ho.at(time)[indiv].end(), locus) != ho.at(time)[indiv].end());
        }

        bool   hidden(size_t time, size_t indiv, size_t locus) const      { return (missing(time, indiv, locus) || hold_out(time, indiv, locus)); }

        bool   hemizygous(size_t time, size_t indiv) const  { return hemi[time][indiv]; }

        bool   has_hold_out() { return ho.size() != 0; }

    private:
        typedef std::vector<std::vector<int> >  matrix;
        const std_vector3<short>                *snps;         // the full SNP data set
        std::map<int, matrix>                   ho;            // 
        std::vector<int>                        sample_gen;    // sample_gen[t] gives the generation number corresponding to time step t.
        std::vector<std::vector<bool> >         hemi;          // hemi[t][d] is true if individual is hemizygous
};

#endif