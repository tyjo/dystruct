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
        SNPData(SNPData& d) : ho(d.ho), sample_gen(d.sample_gen) { this->snps = d.snps; }
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

    private:
        typedef std::vector<std::vector<int> >  matrix;
        const std_vector3<short>                *snps;         // the full SNP data set
        std::map<int, matrix>                   ho;            // 
        std::vector<int>                        sample_gen;    // sample_gen[t] gives the generation number corresponding to time step t.
};

#endif