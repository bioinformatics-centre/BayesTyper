
/*
Sample.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


The MIT License (MIT)

Copyright (c) 2016 Jonas Andreas Sibbesen and Lasse Maretty

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/


#include <fstream>

#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string.hpp"

#include "Sample.hpp"


Sample::Sample(const string & sample_line) {

    vector<string> split_sample_line;
    boost::split(split_sample_line, sample_line, boost::is_any_of("\t"));
    assert(split_sample_line.size() == 3);

    name = split_sample_line.at(0);

    if (split_sample_line.at(1) == "M") {

        gender = Utils::Gender::Male;

    } else {

        assert(split_sample_line.at(1) == "F");
        gender = Utils::Gender::Female;
    }

    file = split_sample_line.at(2);

    ifstream kmer_prefix_infile(file + ".kmc_pre");
    ifstream kmer_suffix_infile(file + ".kmc_suf");

    if (!kmer_prefix_infile.good() or !kmer_suffix_infile.good()) {

        cout << "\nERROR: " << file << ".kmc_pre or " << file << ".kmc_suf does not exist - or you do not have sufficient permissions.\n" << endl;
        exit(1);
    }

    kmer_prefix_infile.close();
    kmer_suffix_infile.close();

    ifstream bloom_data_infile(file + ".bloomData");
    ifstream bloom_meta_infile(file + ".bloomMeta");

    if (!bloom_data_infile.good() or !bloom_meta_infile.good()) {

        cout << "\nERROR: " << file << ".bloomData or " << file << ".bloomMeta does not exist - or you do not have sufficient permissions.\n" << endl;
        exit(1);
    }

    bloom_data_infile.close();
    bloom_meta_infile.close();

}