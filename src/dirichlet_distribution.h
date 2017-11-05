#ifndef DIRICHLET_H
#define DIRICHLET_H

#include <algorithm>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <cassert>
#include <vector>

class dirichlet_distribution
{
    public:
        dirichlet_distribution(const std::vector<double>& alpha) { this->alpha = alpha; }
        
        std::vector<double> operator() (boost::random::mt19937& gen) {
            {
                std::vector<double> sample(alpha.size());
                double s = 0.0;
                for (size_t i = 0; i < sample.size(); ++i) {
                    boost::random::gamma_distribution<double> gamma(this->alpha[i], 1);
                    sample[i] = gamma(gen);
                    s += sample[i];
                }
                assert(s != 0);
                for (size_t i = 0; i < sample.size(); ++i)
                    sample[i] /= s;
                return sample;
            }
        }

    private:
        std::vector<double> alpha;
};

#endif