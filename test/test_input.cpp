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

#include <boost/multi_array.hpp>
#include <iostream>
#include <string>
#include <vector>

using std::cout;
using std::endl;
using std::string;
using std::vector;

#include "../src/snp_data.h"
#include "../src/util.h"
#include "../src/vector_types.h"

void test_input1()
{
    cout << "testing correct input file..." << endl;
    string in_file = "test_input1";
    int nloci = 2;
    std_vector3<short> *snps = new std_vector3<short>;
    vector<int> gen_sampled;
    read_snp_matrix(in_file, snps, gen_sampled, nloci);
    cout << endl;
}


void test_error_input1()
{
    cout << "testing unsorted input file..." << endl;
    string in_file = "error_input1";
    int nloci = 1;
    std_vector3<short> *snps = new std_vector3<short>;
    vector<int> gen_sampled;
    read_snp_matrix(in_file, snps, gen_sampled, nloci);
    cout << endl;
}


void test_error_input2()
{
    cout << "testing input file with bad genotype..." << endl;
    string in_file = "error_input2";
    int nloci = 1;
    std_vector3<short> *snps = new std_vector3<short>;
    vector<int> gen_sampled;
    read_snp_matrix(in_file, snps, gen_sampled, nloci);
    cout << endl;
}


int main(int argc, char* const argv[])
{   
    test_input1();
    //test_error_input1();
    test_error_input2();
    return 0;
}