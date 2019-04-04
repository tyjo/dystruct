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
#include <map>
#include <string>
#include <utility>
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
using std::map;
using std::pair;
using std::setprecision;
using std::string;
using std::vector;

string VERSION = "v1.1.0";

void print_help()
{
    cerr << "Program: dystruct (Dynamic Structure)" << endl;
    cerr << "Version: " << VERSION << endl;
    cerr << "Contact: Tyler Joseph <tjoseph@cs.columbia.edu>" << endl;
    cerr << endl;
    cerr << "Usage:   dystruct [options]" << endl;
    cerr << endl;
    cerr << "Required Arguments:" << endl;
    cerr << "\t-h, --help                  " << "Print this help message." << endl;
    cerr << "\t--input FILE                " << "Path to genotype matrix: a LOCI x INDIVIDUAL matrix of genotypes in the" << endl
         << "                                    EIGENSTRAT genotype format. Each genotype is denoted by either 0, 1, 2, or 9," << endl
         << "                                    where 9 identifies missing entries. There are no spaces between entries." << endl
         << "                                    See https://github.com/DReichLab/EIG/tree/master/CONVERTF for more details" << endl
         << "                                    and converting between standard formats." << endl;
    cerr << "\t--generation-times FILE     " << "Path to generation times corresponding to the input file. A list of generation" << endl
         << "                                    times (one per line) for each sample. Samples are assumed to be in the same" << endl;
    cerr << "                                    order as the columns of the input matrix." << endl; 
    cerr << "\t--output STR                " << "A prefix for output files." << endl;
    cerr << "\t--npops INT                 " << "Number of populations." << endl;
    cerr << "\t--nloci INT                 " << "Number of loci. This should match the number of loci in the input file." << endl;
    cerr << "\t--pop-size INT              " << "Effective population size for all populations." << endl;
    cerr << "\t--seed INT                  " << "Random seed used to initialize variational parameters" << endl;
    cerr << endl;
    cerr << "Optional Arguments:" << endl;
    cerr << "\t--hold-out-fraction DOUBLE  " << "(=0) Optional. Partitions nloci * hold_out_fraction loci into a hold out" << endl
         << "                                    set. The hold out set contains at most one individual per site." << endl;
    cerr << "\t--hold-out-seed INT         " << "(=28149) Optional. Random seed used to partition SNP data into hold out" << endl
         << "                                    and training sets. Use the same seed across replicates to fix the hold" << endl
         << "                                    out set." << endl;
    cerr << "\t--epochs INT                " << "(=50) Optional. Number of epochs to run before terminating." << endl;
    cerr << "\t--no-multi-init             " << "(=true) Optional. Turns off multiple initialization." << endl;
    /*cerr << "\t--labels FILE               " << "Optional. Experimental. Population label file path for supervised analysis." << endl 
         << "                                    Labels should be in {0,...,npops - 1}. One label per line in the same order" << endl
         << "                                    as the input matrix. Individuals without a population assignment should be" << endl
         << "                                    labeled by -1." << endl;*/
}


enum OPTIONS
{
    INPUT,
    GENERATION_TIMES,
    OUTPUT,
    NPOPS,
    NLOCI,
    POP_SIZE,
    SEED,
    HOLD_OUT_FRACTION,
    HOLD_OUT_SEED,
    EPOCHS,
    MULTI_INIT,
    LABELS
};


static struct option long_options[] =
{
    {"input"             , required_argument, NULL, INPUT             },
    {"generation-times"  , required_argument, NULL, GENERATION_TIMES  },
    {"output"            , required_argument, NULL, OUTPUT            },
    {"npops"             , required_argument, NULL, NPOPS             },
    {"nloci"             , required_argument, NULL, NLOCI             },
    {"pop-size"          , required_argument, NULL, POP_SIZE          },
    {"seed"              , required_argument, NULL, SEED              },
    {"hold-out-fraction" , required_argument, NULL, HOLD_OUT_FRACTION },
    {"hold-out-seed"     , required_argument, NULL, HOLD_OUT_SEED     },
    {"epochs"            , required_argument, NULL, EPOCHS            },
    {"no-multi-init"     , no_argument      , NULL, MULTI_INIT        },
    {"labels"            , required_argument, NULL, LABELS            },
    {NULL, no_argument, NULL, 0}
};




int main(int argc, char* const argv[])
{
    if (argc <= 2) {
        print_help();
        return 1;
    }

    string in_file           = "";
    string in_gen_times_file = "";
    string out_file          = "";
    int random_seed          = 0;
    int hold_out_seed        = 28149;
    int npop                 = 0;
    int nloci                = 0;
    double hold_out_fraction = 0;
    double pop_size          = 0;
    int epochs               = 50;
    string label_file        = "";
    bool multi_init          = true;

    int c;
    int option_index;
    while ( (c = getopt_long(argc, argv, "", long_options, &option_index)) != -1 ) {
        switch (c) {
            case INPUT:
                in_file = optarg;
                break;
            case GENERATION_TIMES:
                in_gen_times_file = optarg;
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
            case EPOCHS:
                epochs = atoi(optarg);
                break;
            case MULTI_INIT:
                multi_init = false;
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
    else if (in_gen_times_file == "") {
        cerr << "missing argument: --generation-times" << endl;
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
    else if (epochs < 0) {
        cerr << "--epochs must be greater than 0" << endl;
        return 1;
    }

    // initialize random number generator
    mt19937 gen(random_seed);
    std_vector3<short> *snps = new std_vector3<short>;

    vector<int> gen_sampled;
    map<int, pair<int, int> > sample_map = read_snp_matrix(in_file, in_gen_times_file, snps, gen_sampled, nloci);
    if (hold_out_fraction > 0)
        cout << "constructing hold out set..." << endl;
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

    //cout << "initializing variational parameters..." << endl;
    SVI svi(npop, theta_prior, pop_size, snp_data, gen, nloci, epochs, sample_map, labels, multi_init, use_labels);

    //cout << "running..." << endl;
    svi.run_stochastic();
    svi.write_results(out_file);

    return 0;
}
