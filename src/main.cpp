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
#include <cmath>
#include <getopt.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "svi.h"
#include "snp_data.h"
#include "util.h"

using boost::random::mt19937;
using boost::random::uniform_int_distribution;

using std::cerr;
using std::cout;
using std::exit;
using std::endl;
using std::setprecision;
using std::string;
using std::vector;

string VERSION = "v0.2.4";


void print_help()
{
    cerr << "Program: dystruct (Dynamic Structure)" << endl;
    cerr << "Version: " << VERSION << endl;
    cerr << "Contact: Tyler Joseph <tjoseph@cs.columbia.edu>" << endl;
    cerr << endl;
    cerr << "Usage:   dystruct [options]" << endl;
    cerr << endl;
    cerr << "Required Arguments:" << endl;
    cerr << "\t--input FILE                " << "Genotype file path. An LOCI x INDIVIDUAL matrix of genotypes. The header" << endl
         << "                                    is the sample time in generations. Samples must be ordered in increasing" << endl
         << "                                    generation time. With genotypes encoder by 0, 1, 2, or 9 for missing data." << endl;
    cerr << "\t--output STR                " << "A file prefix for output files." << endl;
    cerr << "\t--npops INT                 " << "Number of populations." << endl;
    cerr << "\t--nloci INT                 " << "Number of loci. This should match the number of loci in the input file." << endl;
    cerr << "\t--pop-size INT              " << "Effective population size for all populations." << endl;
    cerr << "\t--seed INT                  " << "Random seed used to initialize variational parameters" << endl;
    cerr << endl;
    cerr << "Optional Arguments:" << endl;
    cerr << "\t--hold-out-fraction DOUBLE  " << "(=0) Optional. Partitions nloci * hold_out_fraction loci into a hold out" << endl
         << "                                    set. The hold out set contains at most one site per individual." << endl;
    cerr << "\t--hold-out-seed INT         " << "(=28149) Optional. Random seed used to partition SNP data into hold out" << endl
         << "                                    and training sets. Use the same seed across replicates to fix the hold" << endl
         << "                                    out set." << endl;
    cerr << "\t--tol DOUBLE                " << "(=1) Optional. Convergence threshold in number of loci. Programs terminates" << endl
         << "                                    when delta < tol." << endl;
    cerr << "\t--step-size-power DOUBLE    " << "(=-0.6) Optional. Adjusts step size for stochastic variational inference." << endl
         << "                                    step_size = (iteration + 1)^step-size-power. Must be in [-1,-0.5)." << endl;
    cerr << "\t--labels FILE               " << "Optional. Experimental. Population label file path for supervised analysis." << endl 
         << "                                    Labels should be in {0,...,npops - 1}. One label per line in the same order" << endl
         << "                                    as the input matrix. Individuals without a population assignment should be" << endl
         << "                                    labeled by -1." << endl;
}


enum OPTIONS
{
    INPUT,
    OUTPUT,
    NPOPS,
    NLOCI,
    POP_SIZE,
    SEED,
    HOLD_OUT_FRACTION,
    HOLD_OUT_SEED,
    TOL,
    STEP_SIZE_POWER,
    LABELS
};


static struct option long_options[] =
{
    {"input"             , required_argument, NULL, INPUT             },
    {"output"            , required_argument, NULL, OUTPUT            },
    {"npops"             , required_argument, NULL, NPOPS             },
    {"nloci"             , required_argument, NULL, NLOCI             },
    {"pop-size"          , required_argument, NULL, POP_SIZE          },
    {"seed"              , required_argument, NULL, SEED              },
    {"hold-out-fraction" , required_argument, NULL, HOLD_OUT_FRACTION },
    {"hold-out-seed"     , required_argument, NULL, HOLD_OUT_SEED     },
    {"tol"               , required_argument, NULL, TOL               },
    {"step-size-power"   , required_argument, NULL, STEP_SIZE_POWER   },
    {"labels"            , required_argument, NULL, LABELS            },
    {NULL, no_argument, NULL, 0}
};




int main(int argc, char* const argv[])
{

    if (argc <= 1) {
        print_help();
        return 1;
    }

    string in_file           = "";
    string out_file          = "";
    int random_seed          = 0;
    int hold_out_seed        = 28149;
    int npop                 = 0;
    int nloci                = 0;
    double hold_out_fraction = 0;
    double pop_size          = 0;
    double tol               = 1;
    double step_power        = -0.6;
    string label_file        = "";

    int c;
    int option_index;
    while ( (c = getopt_long(argc, argv, "", long_options, &option_index)) != -1 ) {
        switch (c) {
            case INPUT:
                in_file = optarg;
                break;
            case OUTPUT:
                out_file = optarg;
                break;
            case NPOPS:
                npop = atoi(optarg);
                break;
            case NLOCI:
                nloci = atoi(optarg);
                break;
            case POP_SIZE:
                pop_size = atof(optarg);
                break;
            case SEED:
                random_seed = atoi(optarg);
                break;
            case HOLD_OUT_FRACTION:
                hold_out_fraction = atof(optarg);
                break;
            case HOLD_OUT_SEED:
                hold_out_seed = atoi(optarg);
                break;
            case TOL:
                tol = atof(optarg);
                break;
            case STEP_SIZE_POWER:
                step_power = atof(optarg);
                break;
            case LABELS:
                label_file = optarg;
                break;
            default:
                // getopt prints an error
                return 1;
                break;
        }
    }

    // check input
    if (in_file == "") {
        cerr << "missing argument: --input" << endl;
        return 1;
    }
    else if (random_seed == 0) {
        cerr << "missing argument: --seed" << endl;
        return 1;
    }
    else if (random_seed < 0) {
        cerr << "argument error: --seed must be a positive integer" << endl;
        return 1;
    }
    else if (npop == 0) {
        cerr << "missing argument: --npop" << endl;
        return 1;
    }
    else if (nloci == 0) {
        cerr << "missing argument: --nloci" << endl;
        return 1;
    }
    else if (hold_out_fraction < 0 || hold_out_fraction >= 1) {
        cerr << "proportion of held out sites must be between in [0, 1)" << endl;
        return 1;
    }
    else if (tol < 0) {
        cerr << "--tol must be positive" << endl;
        return 1;
    }
    else if (step_power <= -1 || step_power > -0.5) {
        cerr << "power for step size must be in the interval [-1, -0.5)" << endl;
        return 1;
    }
    
    // initialize random number generator
    mt19937 gen(random_seed);
    std_vector3<short> *snps = new std_vector3<short>;

    vector<int> gen_sampled;
    read_snp_matrix(in_file, snps, gen_sampled, nloci);
    SNPData snp_data(snps, gen_sampled, hold_out_fraction, hold_out_seed);
    vector2<int> labels(boost::extents[snp_data.total_time_steps()][snp_data.max_individuals()]);
    bool use_labels = false;
    if (label_file != "") {
        labels = read_pop_labels(label_file, snp_data);
        use_labels = true;
    }
    
    vector<double> theta_prior;
    for (int k = 0; k < npop; ++k) {
       theta_prior.push_back(1.0 / npop);
    }

    cout << "initializing variational parameters..." << endl;
    SVI svi(npop, theta_prior, pop_size, snp_data, gen, nloci, tol, step_power, labels, use_labels);

    cout << "running..." << endl;
    svi.run_stochastic();
    svi.write_results(out_file);
    
    return 0;
}
