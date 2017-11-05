#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <cassert>
#include <map>
#include <set>
#include <vector>
#include "snp_data.h"

using std::pair;
using std::map;
using std::vector;

#include <iostream>

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

    size_t nloci = snps[0][0].size();
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
        int t = tdist(gen);

        // pick an individual
        boost::random::uniform_int_distribution<int> idist(0, snps[t].size() - 1);
        int individual = idist(gen);

        // make sure we have data
        while (missing(t, individual, draw)) {
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