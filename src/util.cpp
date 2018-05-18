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

using std::cerr;
using std::cout;
using std::count;
using std::exit;
using std::endl;
using std::find;
using std::ifstream;
using std::is_sorted;
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

#include "snp_data.h"
#include "util.h"
#include "vector_types.h"


vector<int> read_generations(string fname, vector<int>& gen_sampled)
{
    ifstream input(fname);
    if (!input.is_open()) {
        cerr << "cannot open " << fname << endl;
        exit(1);
    }

    string line;
    getline(input, line);
    vector<int> generations;
    int g;
    istringstream iss(line);
    while (iss >> skipws >> g) {
        generations.push_back(g);
    }

    if (!is_sorted(generations.begin(), generations.end())) {
        cerr << "Input Error (" << fname << "): samples are not ordered by generation time" << endl;
        exit(1);
    }

    // remove duplicate sample times so we can aggregate samples by generation sampled
    insert_iterator<vector<int> > insert_it(gen_sampled, gen_sampled.begin());
    unique_copy(generations.begin(), generations.end(), insert_it);

    input.close();

    return generations;
}


void check_input_file(string fname, int nloci)
{   
    cout << "checking input file..." << endl;
    ifstream input(fname);
    if (!input.is_open()) {
        cerr << "cannot open " << fname << endl;
        exit(1);
    }

    string line;
    int g;
    int n_columns = 0;
    getline(input, line);

    // count the number of columns
    istringstream iss(line);
    while (iss >> skipws >> g) {
        n_columns++;
    }

    int locus_count = 0;
    int col_count = 0;
    while(getline(input, line)) {
        locus_count++;
        iss = istringstream(line);
        col_count = 0;

        while(iss >> skipws >> g) {
            col_count++;

            if (!(g == 0 || g == 1 || g == 2 || g == 9)) {
                cerr << "Input Error (" << fname << "): line " << locus_count << " column "
                     << col_count + 1 << " has an invalid entry." << endl;
                cerr << "Genotypes must be 0, 1, or 2 if known, 9 if missing or unknown." << endl;
                exit(1);
            }
        }

        if (col_count != n_columns) {
            cerr << "Input Error (" << fname << "): line " << locus_count << " has "
                 << col_count << " samples, but header has " << n_columns << "." << endl;
            exit(1);
        }
    }

    if (nloci != locus_count) {
        cerr << "Input Warning (" << fname << "): " << nloci << " were specified, but "
             << locus_count << " were found." << endl;
    }
}


void read_snp_matrix(string fname, std_vector3<short> *snps, vector<int>& gen_sampled, int nloci)
{
    cout << "loading genotype matrix (this should not take more than a few minutes)..." << endl;
    check_input_file(fname, nloci);
    vector<int> generations = read_generations(fname, gen_sampled);

    cout << "found " << generations.size() << " samples at " << gen_sampled.size() << " time points..." << endl;
    cout << "using " << nloci << " loci..." << endl; 

    ifstream input(fname);
    if (!input.is_open()) {
        cerr << "cannot open " << fname << endl;
        exit(1);
    }

    string line;
    getline(input, line); // first line is the header

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
    istringstream iss(line);
    while (getline(input, line)) {
        iss = istringstream(line);
        i = 0;
        map<int, int> gen_to_nsamples; // generation to number of samples
        while (iss >> skipws >> genotype) {
            gen = generations[i];
            index = find(gen_sampled.begin(), gen_sampled.end(), gen) - gen_sampled.begin();
            (*snps)[index][gen_to_nsamples[gen]][l] = genotype;
            gen_to_nsamples[gen]++;
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
    input.close();
    return labels;
}