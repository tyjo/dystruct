#ifndef KALMAN_SMOOTHER_H
#define KALMAN_SMOOTHER_H
#define NDEBUG // for boost

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <vector>
#include <utility>

#include "snp_data.h"
#include "vector_types.h"

namespace ublas = boost::numeric::ublas;

class VariationalKalmanSmoother
{
    public:
        // Takes the current value of the pseudo-outputs, and the variances from the state space model.
        VariationalKalmanSmoother(const SNPData& snp_data,
                                  const vector3<double>& outputs,
                                  double initial_mean,
                                  const vector3<double>& phi,
                                  const vector3<double>& zeta,
                                  double pop_size,
                                  size_t pop,
                                  size_t locus);

        VariationalKalmanSmoother() { }

        // optimizes variational pseudo-outputs at one locus in one population
        // across all time steps using a conjugate gradient algorithm.
        void maximize_pseudo_outputs();

        // computes the marginal mean and marginal variance at one locus in one
        // population across all time steps and sets freqs to these estimates.
        void set_marginals(vector4<double>& freqs, size_t pop, size_t locus);

        void set_outputs(vector3<double>& outputs);

        // computes the forward mean and variance
        void compute_forward_equations();

        // computes the marginal mean and variance using a backward recurrence
        void compute_backward_equations();

        // computes the partial derivatives of the forward means with
        // respect to the pseudo outputs and initial mean
        inline void compute_forward_partials();

        // computes the partial derivatives of the backward means 
        // with respect to the pseudo outputs and initial mean
        inline void compute_backward_partials();

        // compute the terms in the ELBO that depend on the pseudo-outputs
        inline double compute_objective();

        // compute the terms in the ELBO at the given locus that depend on population size
        double get_initial_mean() { return initial_mean; }

    private:
        inline bool converged(const ublas::vector<double>& v1, const ublas::vector<double>& v2);

        ublas::vector<double> f_mean;                         // forward means
        ublas::vector<double> f_var;                          // forward variances
        ublas::vector<double> marg_mean;                      // marginal/backward means
        ublas::vector<double> marg_var;                       // marginal/backward variances
        ublas::vector<double> var;                            // state space variance
        double out_var;                                       // variational parameter for the output variance of the pseudo-outputs
        ublas::matrix<double> forward_partials;
        ublas::matrix<double> partials;                       // partial derivatives of marg_mean[t] with respect to outputs[s] (partials[t][s])
        ublas::vector<double> param_forward_derivative;       // forward derivative with respect to initial mean
        ublas::vector<double> param_backward_derivative;      // derivative with respect to initial mean
        ublas::vector<double> sum_phi;
        ublas::vector<double> sum_zeta;
        double initial_mean;
        double initial_variance;
        double pop_size;
        ublas::vector<double> dt;                             // number of generations between time steps
        size_t total_time_steps;
        ublas::vector<double> parameters;
        ublas::vector<double> phi_zeta;
        ublas::matrix<double> diff_partials;
        ublas::vector<double> diff_means;
        size_t time_points;
        size_t pop;
        size_t locus;
};

#endif