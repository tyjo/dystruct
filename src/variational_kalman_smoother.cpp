#define NDEBUG // for boost
#define BOOST_MATH_PROMOTE_DOUBLE_POLICY false

#include <algorithm>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

#include "snp_data.h"
#include "variational_kalman_smoother.h"

using std::abs;
using std::pair;
using std::cout;
using std::endl;
using std::max;
using boost::numeric::ublas::matrix;
using boost::numeric::ublas::vector;

#include <iomanip>
using std::setprecision;

VariationalKalmanSmoother::VariationalKalmanSmoother(const SNPData& snp_data,
                                                     const vector3<double>& outputs,
                                                     double initial_mean,
                                                     const vector3<double>& phi,
                                                     const vector3<double>& zeta,
                                                     double pop_size,
                                                     size_t pop,
                                                     size_t locus) :
    parameters(snp_data.total_time_steps() + 1, 0),
    var(snp_data.total_time_steps(), 0),
    f_mean(snp_data.total_time_steps(), 0),
    f_var(snp_data.total_time_steps(), 0),
    marg_mean(snp_data.total_time_steps(), 0),
    marg_var(snp_data.total_time_steps(), 0),
    forward_partials(snp_data.total_time_steps(), snp_data.total_time_steps() + 1, 0), // time points by number of parameters
    partials(snp_data.total_time_steps(), snp_data.total_time_steps() + 1, 0),
    param_forward_derivative(snp_data.total_time_steps(), 0),
    param_backward_derivative(snp_data.total_time_steps(), 0),
    sum_phi(snp_data.total_time_steps(), 0),
    sum_zeta(snp_data.total_time_steps(), 0),
    phi_zeta(snp_data.total_time_steps(), 0),
    dt(snp_data.total_time_steps(), 0),
    diff_partials(snp_data.total_time_steps(), snp_data.total_time_steps()+ 1, 0),
    diff_means(snp_data.total_time_steps(), 0),
    time_points(snp_data.total_time_steps()),
    pop(pop),
    locus(locus)

{
    this->initial_mean = initial_mean;
    this->pop_size = pop_size;
    this->out_var = 0.001;
    this->total_time_steps = snp_data.total_time_steps();
    this->dt[0] = 1;
    this->var[0] = 1. / (12*pop_size);
    this->parameters[time_points] = initial_mean;
    this->initial_variance = var[0];
    for (size_t t = 0; t < time_points; ++t) {
        parameters[t] = outputs[pop][locus][t];

        if (t > 0) {
            // store generation time
            this->dt[t] = snp_data.get_sample_gen(t) - snp_data.get_sample_gen(t-1);
        
            // compute state space variance
            var[t] = dt[t] / (12.*pop_size);
        }
        
        // sums used in gradient/objective function calculation
        // stored here so we don't have to recompute every iteration
        // during optimization
        sum_phi[t] = 0;
        sum_zeta[t] = 0;
        for (size_t d = 0; d < snp_data.total_individuals(t); ++d) {
            if (snp_data.hidden(t, d, locus)) continue;
            sum_phi[t] += snp_data.genotype(t, d, locus) * phi[t][d][pop];
            sum_zeta[t] += (2 - snp_data.genotype(t, d, locus)) * zeta[t][d][pop];
        }
    }
}



void VariationalKalmanSmoother::maximize_pseudo_outputs()
{
    // store the previous and next iteration's parameter values
    // used to check for convergence. The last element of the vector
    // is the initial mean.
    vector<double> prev(time_points + 1, 0); 
    vector<double> next(time_points + 1, 0);

    // vector of partial derivatives of marginal means with respect to pseudo-outputs
    // the last component is the partial derivative with respect to the initial mean
    vector<double> grad(time_points + 1, 0);

    // the vectors p are the conjugate vectors
    vector<double> p(time_points + 1, 0);

    double b = 0;
    double prev_norm = 0;
    double norm = 0;
    double step_size = 0.0001;

    // the current and future values of the objective function
    double obj = 0.0;
    double nxt = 0.0;

    unsigned int it = 0;
    unsigned int max_iter = max(2*(time_points + 1), 10ul);
    //unsigned int max_iter = 50000;

    double min_out = 0.01;
    double max_out = 1 - min_out;

    while ( (!converged(prev, next)  || it < 2) && it < max_iter) {
        compute_forward_partials();
        compute_backward_partials();
        it++;

        grad = prod(diff_means, diff_partials) + prod(phi_zeta, partials);
        prev_norm = norm;
        norm = norm_2(grad);

        //
        // conjugate gradient update
        //

        if (norm < 0.01 && norm != 0)
            break;

        if (norm != 0 && prev_norm != 0)
            b = norm / prev_norm;

        p = grad + b*p;
        prev = parameters;

        // Wolfe conditions
        // we are minimizing the negative of the objective
        step_size = 0.0001;
        obj = compute_objective();
        nxt = 0;
        while ((-nxt > -obj - 0.0001*step_size*inner_prod(grad, p)) || nxt == 0) {
            if (nxt != 0)
                step_size = step_size / 2; 
            for (size_t s = 0; s < parameters.size(); ++s) {
                parameters[s] = prev[s] + step_size*p[s];
                if (parameters[s] >= max_out) parameters[s] = max_out;
                if (parameters[s] <= min_out) parameters[s] = min_out;
            }
            initial_mean = parameters[time_points];
            nxt = compute_objective();
        }
        
        next = parameters;
    }
}



bool VariationalKalmanSmoother::converged(const vector<double>& v1, const vector<double>& v2)
{
    bool has_converged = true;
    for (size_t i = 0; i < v1.size(); ++i) {
        if (abs(v1[i] - v2[i]) > 10E-10) {
            has_converged = false;
            break;
        }
    }
    return has_converged;
}



void VariationalKalmanSmoother::set_marginals(vector4<double>& freqs, size_t k, size_t l)
{
    compute_forward_equations();
    compute_backward_equations();
    for (size_t t = 0; t < time_points; ++t) {
        freqs[t][k][l][0] = marg_mean[t];
        freqs[t][k][l][1] = marg_var[t];
    }
}


void VariationalKalmanSmoother::set_outputs(vector3<double>& outputs)
{
    for (size_t t = 0; t < time_points; ++t) {
        outputs[pop][locus][t] = parameters[t];
    }
}



void VariationalKalmanSmoother::compute_forward_equations()
{
    f_mean[0] = (out_var / (initial_variance + var[0] + out_var))*initial_mean
              + (1 - out_var / (initial_variance + var[0] + out_var))*parameters[0];
    f_var[0] = (out_var / (initial_variance + var[0] + out_var)) * (initial_variance + var[0]);
    for (size_t t = 1; t < time_points; ++t) {
        f_mean[t] = (out_var / (f_var[t-1] + var[t] + out_var))*f_mean[t-1]
                  + (1 - out_var / (f_var[t-1] + var[t] + out_var))*parameters[t];
        f_var[t] = (out_var / (f_var[t-1] + var[t] + out_var)) * (f_var[t-1] + var[t]);
    }
}



void VariationalKalmanSmoother::compute_backward_equations()
{
    marg_mean[time_points - 1] = f_mean[time_points - 1];
    marg_var[time_points - 1] = f_var[time_points - 1];
    for (size_t t = time_points - 2; t < time_points; --t) {
        marg_mean[t] = (var[t] / (f_var[t] + var[t]))*f_mean[t]
                     + (1 - var[t]/ (f_var[t] + var[t]))*marg_mean[t+1];
        marg_var[t] = f_var[t] + (f_var[t] / (f_var[t] + var[t]))*(f_var[t] / (f_var[t] + var[t])) * (marg_var[t+1] - f_var[t] - var[t]);
    }
}



void VariationalKalmanSmoother::compute_forward_partials()
{   
    // forward means
    f_mean[0] = (out_var / (initial_variance + var[0] + out_var))*initial_mean
              + (1 - out_var / (initial_variance + var[0] + out_var))*parameters[0];
    f_var[0] = (out_var / (initial_variance + var[0] + out_var)) * (initial_variance + var[0]);
    
    // initial mean
    forward_partials(0, time_points) = (out_var / (initial_variance + var[0] + out_var));

    forward_partials(0, 0) = (1 - out_var / (initial_variance + var[0] + out_var));
    for (size_t t = 1; t < time_points; ++t) {

        // forward means
        f_mean[t] = (out_var / (f_var[t-1] + var[t] + out_var))*f_mean[t-1]
                  + (1 - out_var / (f_var[t-1] + var[t] + out_var))*parameters[t];
        f_var[t] = (out_var / (f_var[t-1] + var[t] + out_var)) * (f_var[t-1] + var[t]);
        
        for (size_t s = 0; s < time_points; ++s) {
             forward_partials(t, s) = (out_var / (f_var[t-1] + var[t] + out_var)) * forward_partials(t-1, s);
             if (s == t)
                forward_partials(t, s) += (1 - out_var/(f_var[t-1] + var[t] + out_var));
        }
        // initial mean
        forward_partials(t, time_points) = forward_partials(t-1, time_points) * (out_var / (f_var[t-1] + var[t] + out_var));
    }
}



void VariationalKalmanSmoother::compute_backward_partials()
{   
    // marginal means
    marg_mean[time_points - 1] = f_mean[time_points - 1];
    marg_var[time_points - 1] = f_var[time_points - 1];

    // initial mean
    partials(time_points - 1, time_points) = forward_partials(time_points - 1, time_points);
    diff_partials(time_points - 1, time_points) = partials(time_points - 1, time_points);
    for (size_t t = time_points - 2; t < time_points; --t) {
        // marginal means
        marg_mean[t] = (var[t] / (f_var[t] + var[t]))*f_mean[t]
                     + (1 - var[t]/ (f_var[t] + var[t]))*marg_mean[t+1];
        marg_var[t] = f_var[t] + pow(f_var[t] / (f_var[t] + var[t]), 2) * (marg_var[t+1] - f_var[t] - var[t]);
        
        partials(t, time_points) = forward_partials(t, time_points) * var[t] / (f_var[t] + var[t]) 
                                   + partials(t+1, time_points) * (1 - var[t] / (f_var[t] + var[t]));
        diff_partials(t, time_points) = partials(t, time_points);
        diff_partials(t+1, time_points) -= partials(t, time_points);
    }
    diff_partials(0, time_points) -= 1;

    for (size_t s = 0; s < time_points; ++s) {
        partials(time_points - 1, s) = forward_partials(time_points - 1, s);
        diff_partials(time_points - 1, s) = partials(time_points - 1, s);
        for (size_t t = time_points - 2; t < time_points; --t) {
            partials(t, s) = (var[t] / (f_var[t] + var[t]))*forward_partials(t, s)
                          + (1 - var[t] / (f_var[t] + var[t]))*partials(t+1, s);

            diff_partials(t, s) = partials(t, s);
            diff_partials(t+1, s) -= partials(t, s);
        }

        // precompute difference of means and difference of partials
        diff_means[s] = marg_mean[s];
        if (s > 0) {
            diff_means[s] -= marg_mean[s-1];
        }
        else {
            diff_means[s] -= initial_mean;
        }
        diff_means[s] *= -12*pop_size/dt[s];
        phi_zeta[s] = sum_phi[s]*(1./marg_mean[s] + marg_var[s]/(marg_mean[s]*marg_mean[s]*marg_mean[s])) + sum_zeta[s] / (marg_mean[s] - 1);
    }
}



double VariationalKalmanSmoother::compute_objective()
{
    compute_forward_equations();
    compute_backward_equations();

    double obj = 0.0;
    double m0 = 0;
    double m1 = 0;
    for (size_t t = 0; t < time_points; ++t) {
        if (t > 0) {
            m0 = marg_mean[t-1];
        }
        else {
            m0 = initial_mean;
        }
        m1 = marg_mean[t];
        // E[log p(beta)] - E[log q(beta)]
        obj += -6*pop_size*(m1 - m0)*(m1 - m0)/dt[t] - log(dt[t] / (12*pop_size)) - 12*pop_size*marg_var[t]/dt[t];
        obj += 0.5*log(marg_var[t]);

        // E[log p(x)]
        obj += sum_phi[t]*(log(m1) - marg_var[t]/(2*m1*m1)) + sum_zeta[t]*log(1 - m1);
    }
    obj += 6*pop_size*marg_var[time_points - 1]/dt[time_points - 1];
    return obj;
}
