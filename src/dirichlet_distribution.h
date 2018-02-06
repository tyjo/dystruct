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