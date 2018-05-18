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

#define NDEBUG
#define BOOST_MATH_PROMOTE_DOUBLE_POLICY false

#include <algorithm>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <cassert>
#include <cmath>
#include <random>
#include <string>
#include <ios>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <utility>

#include "cavi.h"
#include "dirichlet_distribution.h"
#include "snp_data.h"
#include "variational_kalman_smoother.h"
#include "vector_types.h"

#include <sstream>
using std::ifstream;
using std::istringstream;
using std::setw;
using std::skipws;

using boost::math::digamma;
using boost::math::lgamma;
using boost::random::mt19937;
using boost::random::gamma_distribution;
using boost::random::normal_distribution;
using boost::random::uniform_real_distribution;
using boost::random::uniform_int_distribution;

using std::abs;
using std::cout;
using std::cerr;
using std::endl;
using std::flush;
using std::ios;
using std::max;
using std::min;
using std::ofstream;
using std::pair;
using std::setprecision;
using std::setfill;
using std::string;
using std::vector;


Cavi::Cavi(int                     npops,
           std::vector<double>     mixture_prior,
           double                  pop_size,
           SNPData                 snp_data,
           boost::random::mt19937& gen,
           size_t                  nloci,
           double                  tol,
           double                  step_power,
           vector2<int>            labels,
           bool                    using_labels) :
           npops(npops),
           nloci(nloci),
           nsteps(snp_data.total_time_steps()),
           initial_freq(boost::extents[npops][nloci]),
           freqs(boost::extents[snp_data.total_time_steps()][npops][nloci][2]),
           pseudo_outputs(boost::extents[npops][nloci][snp_data.total_time_steps()]),
           phi(boost::extents[snp_data.total_time_steps()][snp_data.max_individuals()][npops]),
           zeta(boost::extents[snp_data.total_time_steps()][snp_data.max_individuals()][npops]),
           theta(boost::extents[snp_data.total_time_steps()][snp_data.max_individuals()][npops]),
           nloci_indv(boost::extents[snp_data.total_time_steps()][snp_data.max_individuals()]),
           sample_iter(boost::extents[snp_data.total_time_steps()][snp_data.max_individuals()]),
           tol(tol),
           step_power(step_power),
           labels(labels)
{   
    this->mixture_prior = mixture_prior;
    this->pop_size = pop_size;
    this->snp_data = snp_data;
    this->gen = gen;
    this->using_labels = using_labels;

    cout << setprecision(5);
    initialize_variational_parameters();
}



void Cavi::initialize_variational_parameters()
{
    for (size_t t = 0; t < nsteps; ++t) {
        for (size_t d = 0; d < snp_data.total_individuals(t); ++d) {
            nloci_indv[t][d] = 0;
            for (size_t l = 0; l < nloci; ++l) {
                if (!snp_data.hidden(t, d, l))
                    nloci_indv[t][d] += 1;
                sample_iter[t][d] = 0;
            }
        }
    }

    // initial allele frequency
    uniform_real_distribution<double> uniform(0.2, 0.8);
    for (size_t k = 0; k < npops; ++k) {
        for (size_t l = 0; l < nloci; ++l) {
            initial_freq[k][l] = uniform(gen);
        }
    }

    for (size_t t = 0; t < nsteps; ++t) {
        for (size_t k = 0; k < npops; ++k) {
            for (size_t l = 0; l < nloci; ++l) {
                freqs[t][k][l][0] = initial_freq[k][l];
                freqs[t][k][l][1] =  1.0/(12*pop_size);
            }
        }

        for (size_t d = 0; d < snp_data.total_individuals(t); ++d) {
            gamma_distribution<double> gamma(0.65*(double)nloci_indv[t][d]/100, 100);
            for (size_t k = 0; k < npops; ++k) {
                if (using_labels && labels[t][d] != -1 && labels[t][d] != (int)k) {
                    theta[t][d][k] = 1;
                }
                else {
                    theta[t][d][k] = gamma(gen);
                }


                phi[t][d][k] = 1.0 / npops;
                zeta[t][d][k] = 1.0 / npops;
            }
        }
    }

    // initialize pseudo-outputs
    for (size_t k = 0; k < npops; ++k) {
        for (size_t l = 0; l < nloci; ++l) {
            for (size_t t = 0; t < nsteps; ++t) {
                pseudo_outputs[k][l][t] = uniform(gen);
            }
        }
    }
}



double Cavi::compute_ho_log_likelihood()
{
    double log_lk = 0;
    double p = 0;
    double s = 0;

    for (size_t t = 0; t < nsteps; ++t) {
        for (size_t d = 0; d < snp_data.total_individuals(t); ++d) {
            for (size_t l = 0; l < nloci; ++l) {
                if (!snp_data.hold_out(t,d,l) || snp_data.missing(t,d,l)) continue;
                
                load_auxiliary_parameters(l);
                update_allele_frequencies(l);

                s = 0;
                for (size_t k = 0; k < npops; ++k)
                    s += theta[t][d][k];
                
                p = 0;
                for (size_t k = 0; k < npops; ++k)
                    p += (theta[t][d][k]/s)*freqs[t][k][l][0];

                if (p >= 1) p = 0.999;
                else if (p == 0) p = 0.001;

                double m = max(snp_data.genotype(t,d,l), 2 - snp_data.genotype(t,d,l));
                log_lk += log(2.0) -  log(m) + snp_data.genotype(t,d,l)*log(p) + (2 - snp_data.genotype(t,d,l))*log(1-p);
            }
        }
    }
    return log_lk;
}



bool Cavi::update_auxiliary_parameters(int sample)
{
    bool converged = true;
    for (size_t t = 0; t < nsteps; ++t) {
        
        // better to parallelize inner loop since sample
        // tend to be sparse in t but dense in d
        #pragma omp parallel for reduction(&&:converged)
        for (size_t d = 0; d < snp_data.total_individuals(t); ++d) {
            if (snp_data.hidden(t,d,sample)) continue;
            vector<double> prev_phi(npops);
            vector<double> prev_zeta(npops);
            
            for (size_t k = 0; k < npops; ++k) {
                prev_phi[k] = phi[t][d][k];
                prev_zeta[k] = zeta[t][d][k];
            }
            update_auxiliary_local(t,d,sample);
            
            for (size_t k = 0; k < npops; ++k) {
                converged = (abs(prev_zeta[k] - zeta[t][d][k]) < 0.01) && (abs(prev_phi[k] - phi[t][d][k]) < 0.01) && converged;
            }
        }
    }
    return converged;
}



void Cavi::load_auxiliary_parameters(int sample)
{
    for (size_t t = 0; t < nsteps; ++t) {
        for (size_t d = 0; d < snp_data.total_individuals(t); ++d) {
            if (snp_data.hidden(t,d,sample)) continue;
            update_auxiliary_local(t,d,sample);
        }
    }
}



void Cavi::update_auxiliary_local(size_t t, size_t d, size_t l)
{
    double dgma = 0.0;
    double max_phi = phi[t][d][0];
    double max_zeta = zeta[t][d][0];


    // computes log of auxiliary parameters (log-sum-exp trick for underflow)
    for (size_t k = 0; k < npops; ++k) {
        if (using_labels && labels[t][d] != -1 && labels[t][d] != (int)k) {
            phi[t][d][k] = log(1E-8);
            zeta[t][d][k] = log(1E-8);
        }
        else if (using_labels && labels[t][d] == (int)k) {
            phi[t][d][k] = log(1 - (npops - 1)*1E-4);
            zeta[t][d][k] = log(1 - (npops - 1)*1E-4);
        }
        else {
            dgma = digamma(theta[t][d][k]);     
            phi[t][d][k] = dgma + log(freqs[t][k][l][0]) - freqs[t][k][l][1] / (2*freqs[t][k][l][0]*freqs[t][k][l][0]);
            zeta[t][d][k] = dgma + log(1 - freqs[t][k][l][0]);
        }

        // needed for normalizing constants
        if (phi[t][d][k] > max_phi)
            max_phi = phi[t][d][k];
        if (zeta[t][d][k] > max_zeta)
            max_zeta = zeta[t][d][k];
    }

    // normalizing constants
    // logC = - log \sum_k exp phi[t][d][l][k]
    double sum_exp_phi = 0.0;
    double sum_exp_zeta = 0.0;
    for (size_t k = 0; k < npops; ++k) {
        sum_exp_phi += exp(phi[t][d][k] - max_phi);
        sum_exp_zeta += exp(zeta[t][d][k] - max_zeta);
    }
    double logC_phi = max_phi + log(sum_exp_phi);
    double logC_zeta = max_zeta + log(sum_exp_zeta);

    // normalize
    for (size_t k = 0; k < npops; ++k) {
        phi[t][d][k] = exp(phi[t][d][k] - logC_phi);
        zeta[t][d][k] = exp(zeta[t][d][k] - logC_zeta);

        if (phi[t][d][k] == 0) phi[t][d][k] = 1E-4;
        if (zeta[t][d][k] == 0) zeta[t][d][k] = 1E-4;
    }
}



void Cavi::update_allele_frequencies(int locus)
{
    vector<VariationalKalmanSmoother> vks;
    for (size_t k = 0; k < npops; ++k) {
        vks.push_back(VariationalKalmanSmoother(snp_data, pseudo_outputs, initial_freq[k][locus], phi, zeta, pop_size, k, locus));
    }
    
    #pragma omp parallel for
    for (size_t k = 0; k < npops; ++k) {
        vks[k].maximize_pseudo_outputs();
    }

    for (size_t k = 0; k < npops; ++k) {
        vks[k].set_marginals(freqs, k, locus);
        vks[k].set_outputs(pseudo_outputs);
        initial_freq[k][locus] = vks[k].get_initial_mean();
    }
}



void Cavi::update_mixture_proportions(int locus)
{
    int l = locus;
    double step_size;
    for (size_t t = 0; t < nsteps; ++ t) {
        for (size_t d = 0; d < snp_data.total_individuals(t); ++d) {
            for (size_t k = 0; k < npops; ++k) {
                double update = 0;
                if (snp_data.hidden(t, d, l)) continue;

                step_size = pow(sample_iter[t][d] + 1, step_power);
                
                update += snp_data.genotype(t, d, l)*phi[t][d][k] 
                                + (2 - snp_data.genotype(t, d, l))*zeta[t][d][k];
                theta[t][d][k] += step_size * (mixture_prior[k] + 
                                               nloci_indv[t][d] * update -
                                               theta[t][d][k]
                                              );
                theta[t][d][k] = max(theta[t][d][k], 1.0);
            }
            if (!snp_data.hidden(t, d, l))
                sample_iter[t][d]++;
        }
    }
}


void Cavi::run_stochastic()
{
    int          locus = 0;
    unsigned int it    = 0;
    double ss          = 1;
    bool   converged   = false;

    // monitor convergence
    vector3<double> prev_theta = theta;
    pair<bool, double> theta_converged(false, 100);
    uniform_int_distribution<int> idist(0, nloci - 1);

    while (it < nloci || !theta_converged.first) {
        it++;
        locus = idist(gen);
        ss = pow(it + 1, step_power);

        converged = false;
        while (!converged) {
            converged = update_auxiliary_parameters(locus);
            update_allele_frequencies(locus);
        }
        update_mixture_proportions(locus);

        if (it % 1000 == 0 || it == 1) {
            theta_converged = check_theta_convergence(prev_theta);
            prev_theta = theta;
            cout << "it:\t" << it << setw(10);
            cout << "\tstep size:\t" << ss << setw(10);
            cout << "\tdelta:\t" << theta_converged.second << endl;
            write_temp();
        }
    }
    cout << "hold out log likelihood:\t" << compute_ho_log_likelihood() << endl;
}




pair<bool, double> Cavi::check_theta_convergence(const vector3<double>& prev_theta)
{
    double delta = 0;
    double count = 0;
    for (size_t t = 0; t < nsteps; ++t) {
        for (size_t d = 0; d < snp_data.total_individuals(t); ++d) {
            if (using_labels && labels[t][d] != -1) continue;
            
            for (size_t k = 0; k < npops; ++k) {
                delta += abs(prev_theta[t][d][k] - theta[t][d][k]);
                count += 1;
            }
        }
    }
    delta /= count;
    return pair<bool, double>(delta < tol, delta);
}



void Cavi::write_results(string out_file)
{
    ofstream out_freq(out_file + "_freqs");
    for (size_t t = 0; t < nsteps; ++t) {
        out_freq << t << endl;
        for (size_t l = 0; l < nloci; ++l) {
            for (size_t k = 0; k < npops - 1; ++k) {
                out_freq << freqs[t][k][l][0] << "\t";
            }
            out_freq << freqs[t][npops-1][l][0] << endl;
        }
    }
    out_freq.close();

    ofstream out_theta(out_file + "_theta");
    for (size_t t = 0; t < nsteps; ++t) {
        for (size_t d = 0; d < snp_data.total_individuals(t); ++d) {
            double s = 0;
            for (size_t k = 0; k < npops; ++k) {
                s += theta[t][d][k];
            }
            for (size_t k = 0; k < npops - 1; ++k) {
                out_theta << theta[t][d][k] / s << "\t";
            }
            out_theta << theta[t][d][npops - 1] / s << endl;
        }
    }
    out_theta.close();
}



void Cavi::write_temp()
{
    ofstream out_theta("temp_theta");
    for (size_t t = 0; t < nsteps; ++t) {
        for (size_t d = 0; d < snp_data.total_individuals(t); ++d) {
            for (size_t k = 0; k < npops - 1; ++k) {
                out_theta << theta[t][d][k] << "\t";
            }
            out_theta << theta[t][d][npops - 1] << endl;
        }
    }
    out_theta.close();
}