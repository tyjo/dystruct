#include <boost/program_options.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cmath>
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
using std::setw;
using std::skipws;
using std::string;
using std::unique_copy;
using std::vector;


void read_snp_matrix(string fname, std_vector3<short> *snps, vector<int>& gen_sampled, int nloci);
vector2<int> read_pop_labels(string fname, SNPData& snp_data);

int main(int argc, char* const argv[])
{
    namespace po = boost::program_options;
    po::options_description dystruct("Dystruct Program Usage");
    dystruct.add_options()
        ("help"               , "Print help message")
        ("input,i"            , po::value<string>()         ->required(), "Genotype file path. An L x D matrix of genotypes, where the header is the sample time in generations." 
                                                                          " Samples must be ordered in increasing generation time.")
        ("output,o"           , po::value<string>()         ->default_value(""), "A file prefix that is appended to all output files.")
        ("npops,k"            , po::value<int>()            ->required(), "Number of populations.")
        ("nloci,l"            , po::value<int>()            ->required(), "Number of loci.")
        ("pop_size,z"         , po::value<double>()         ->required()/*->default_value(0)*/, "Specifies population size for all populations.")
        ("seed,s"             , po::value<int>()            ->required(), "Random seed used to initialize variational parameters")
        ("hold_out_fraction,f", po::value<double>()         ->default_value(0), "Optional. Fraction of loci to hold out.")
        ("hold_out_seed,h"    , po::value<int>()            ->default_value(28149), "Optional. Random seed used to partition SNP data into hold out and training sets. Use the same seed across replicates to keep the hold out set fixed.")
        ("labels,b"           , po::value<string>()         ->default_value(""), "Optional. Population labels for supervised analysis. Labels should be in {0,...,npops - 1}. One label per line in the same order as the input matrix. Individuals without a population assignment should be labeled by -1.");

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
    string label_file        = vm["labels"].as<string>();

    // check input
    if (hold_out_fraction < 0 || hold_out_fraction >= 1) {
        cerr << "proportion of held out sites must be between in [0, 1)" << endl;
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
    Cavi cavi(npop, theta_prior, pop_size, snp_data, gen, nloci, labels);

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