#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <sstream>
#include <set>
#include <iomanip>
#include <map>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "cavi.h"
#include "snp_data.h"

using boost::random::mt19937;
using boost::random::uniform_int_distribution;

using std::cerr;
using std::cout;
using std::count;
using std::exit;
using std::endl;
using std::find;
using std::ifstream;
using std::istringstream;
using std::insert_iterator;
using std::map;
using std::pair;
using std::set;
using std::setprecision;
using std::setw;
using std::skipws;
using std::string;
using std::unique_copy;
using std::vector;


void read_snp_matrix(string fname, std_vector3<short> *snps, vector<int>& gen_sampled, int nloci);
vector2<int> read_pop_labels(string fname, SNPData& snp_data);
void print_help();

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
            case STEP_SIZE_POWER:
                step_power = atof(optarg);
                break;
            case LABELS:
                label_file = optarg;
                break;
            default:
                return 1;
                break;
        }
    }

    /*namespace po = boost::program_options;
    po::options_description dystruct("Dystruct Program Usage");
    dystruct.add_options()
        ("help"               , "Print help message")
        ("input"            , po::value<string>()         ->required(), "Genotype file path. An LOCI x INDIVIDUAL matrix of genotypes; the header is the sample time in generations." 
                                                                          " Samples must be ordered in increasing generation time.")
        ("output"           , po::value<string>()         ->default_value(""), "A file prefix for output files.")
        ("npops"            , po::value<int>()            ->required(), "Number of populations.")
        ("nloci"            , po::value<int>()            ->required(), "Number of loci. This should match the number of loci in the input file.")
        ("pop_size"         , po::value<double>()         ->required(), "Specifies population size for all populations.")
        ("seed"             , po::value<int>()            ->required(), "Random seed used to initialize variational parameters")
        ("hold_out_fraction", po::value<double>()         ->default_value(0), "Optional. Partitions nloci * hold_out_fraction loci into a hold out set. The hold out set contains at most one site per individual.")
        ("hold_out_seed"    , po::value<int>()            ->default_value(28149), "Optional. Random seed used to partition SNP data into hold out and training sets. Use the same seed across replicates to fix the hold out set.")
        ("step_size_power"  , po::value<double>()         ->default_value(-0.6, "-0.6"), "Optional. Adjusts step size for stochastic variational inference. step_size = (iteration - offset)^step_power after the first 10000 iterations. The offset ensures the step size does jump between iteration 10000 and 10001. Must be in [-1,-0.5).")
        ("labels"           , po::value<string>()         ->default_value(""), "Optional. Experimental. Population label file path. Population labels for supervised analysis. Labels should be in {0,...,npops - 1}. One label per line in the same order as the input matrix. Individuals without a population assignment should be labeled by -1.");

    if (argc <= 2) {
        cout << dystruct << endl;
        return 1;
    }
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, dystruct), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << dystruct << endl;
        return 1;
    }
    
    string in_file           = vm["input"].as<string>();
    string out_file          = vm["output"].as<string>();
    int random_seed          = vm["seed"].as<int>();
    int hold_out_seed        = vm["hold_out_seed"].as<int>();
    int npop                 = vm["npops"].as<int>(); 
    int nloci                = vm["nloci"].as<int>(); 
    double hold_out_fraction = vm["hold_out_fraction"].as<double>();
    double pop_size          = vm["pop_size"].as<double>();
    double step_power        = vm["step_size_power"].as<double>();
    string label_file        = vm["labels"].as<string>();*/

    // check input
    if (hold_out_fraction < 0 || hold_out_fraction >= 1) {
        cerr << "proportion of held out sites must be between in [0, 1)" << endl;
        return 1;
    }

    if (step_power < -1 || step_power >= -0.5) {
        cerr << "power for step size must be in the interval [-1, -0.5)" << endl;
        return 1;
    }
    
    // initialize random number generator
    cout << "loading genotype matrix..." << endl;
    mt19937 gen(random_seed);
    std_vector3<short> *snps = new std_vector3<short>;

    vector<int> gen_sampled;
    read_snp_matrix(in_file, snps, gen_sampled, nloci);
    SNPData snp_data(snps, gen_sampled, hold_out_fraction, hold_out_seed);
    vector2<int> labels(boost::extents[snp_data.total_time_steps()][snp_data.max_individuals()]);
    if (label_file != "")
        labels = read_pop_labels(label_file, snp_data);
    
    vector<double> theta_prior;
    for (int k = 0; k < npop; ++k) {
       theta_prior.push_back(1.0 / npop);
    }

    cout << "initializing variational parameters..." << endl;
    Cavi cavi(npop, theta_prior, pop_size, snp_data, gen, nloci, step_power, labels);

    cout << "running..." << endl;
    cavi.run_stochastic();
    cavi.write_results(out_file);
    
    return 0;
}


void read_snp_matrix(string fname, std_vector3<short> *snps, vector<int>& gen_sampled, int nloci)
{
    ifstream input(fname);
    if (!input.is_open()) {
        cerr << "cannot find " << fname << endl;
        exit(1);
    }

    string line;

    // read the first line, which is the generation sampled
    getline(input, line);
    vector<int> generations;
    int g;
    istringstream iss(line);
    while (iss >> skipws >> g) {
        generations.push_back(g);
    }
    insert_iterator<vector<int> > insert_it(gen_sampled, gen_sampled.begin());

    // remove duplicate sample times so we can aggregate samples by generation sampled
    unique_copy(generations.begin(), generations.end(), insert_it);

    int ct;
    for (size_t t = 0; t < gen_sampled.size(); ++t) {
        ct = count(generations.begin(), generations.end(), gen_sampled[t]);
        (*snps).push_back(vector2<short>(boost::extents[ct][nloci]));
    }

    short genotype;
    int i;
    int l = 0;
    int gen;
    int index;
    while (getline(input, line)) {
        iss = istringstream(line);
        i = 0;
        map<int, int> gen_samples; // generation to number of samples
        while (iss >> skipws >> genotype) {
            gen = generations[i];
            index = find(gen_sampled.begin(), gen_sampled.end(), gen) - gen_sampled.begin();
            (*snps)[index][gen_samples[gen]][l] = genotype;
            gen_samples[gen]++;
            i++;
        }
        l++;
        if (l == nloci)
            break;
    }
    input.close();
}



vector2<int> read_pop_labels(string fname, SNPData& snp_data)
{
    vector2<int> labels(boost::extents[snp_data.total_time_steps()][snp_data.max_individuals()]);

    ifstream input(fname);
    if (!input.is_open()) {
        cerr << "cannot open label file " << fname << endl;
        exit(1);
    }

    string line;
    istringstream iss;
    int label;
    for (size_t t = 0; t < snp_data.total_time_steps(); ++t) {
        for (size_t d = 0; d < snp_data.total_individuals(t); ++d) {
            getline(input, line);
            iss = istringstream(line);
            iss >> skipws >> label;
            labels[t][d] = label;
        }
    }
    return labels;
}



void print_help()
{
    cerr << "Program: dystruct (Dynamic Structure)" << endl;
    cerr << "Contact: Tyler Joseph <tjoseph@cs.columbia.edu>" << endl;
    cerr << endl;
    cerr << "Usage:   dystruct [options]" << endl;
    cerr << endl;
    cerr << "Options:" << endl;
    cerr << "\t--input FILE                " << "Genotype file path. An LOCI x INDIVIDUAL matrix of genotypes. The header" << endl
         << "                                    is the sample time in generations. Samples must be ordered in increasing" << endl
         << "                                    generation time." << endl;
    cerr << "\t--output STR                " << "A file prefix for output files." << endl;
    cerr << "\t--npops INT                 " << "Number of populations." << endl;
    cerr << "\t--loci INT                  " << "Number of loci. This should match the number of loci in the input file." << endl;
    cerr << "\t--pop-size INT              " << "Effective population size for all populations." << endl;
    cerr << "\t--seed INT                  " << "Random seed used to initialize variational parameters" << endl;
    cerr << "\t--hold-out-fraction DOUBLE  " << "(=0) Optional. Partitions nloci * hold_out_fraction loci into a hold out" << endl
         << "                                    set. The hold out set contains at most one site per individual." << endl;
    cerr << "\t--hold-out-seed INT         " << "(=28149) Optional. Random seed used to partition SNP data into hold out" << endl
         << "                                    and training sets. Use the same seed across replicates to fix the hold" << endl
         << "                                    out set." << endl;
    cerr << "\t--step-size-power DOUBLE    " << "(=-0.6) Optional. Adjusts step size for stochastic variational inference." << endl
         << "                                    step_size = (iteration - offset)^step_power after the first 10000" << endl
         << "                                    iterations. The offset ensures the step size does jump between iteration" << endl
         << "                                    10000 and 10001. Must be in [-1,-0.5)." << endl;
    cerr << "\t--labels FILE               " << "Optional. Experimental. Population label file path for supervised analysis." << endl 
         << "                                    Labels should be in {0,...,npops - 1}. One label per line in the same order" << endl
         << "                                    as the input matrix. Individuals without a population assignment should be" << endl
         << "                                    labeled by -1." << endl;
}