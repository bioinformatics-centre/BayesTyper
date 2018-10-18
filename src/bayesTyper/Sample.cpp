
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

    vector<string> sample_line_split;
    boost::split(sample_line_split, sample_line, boost::is_any_of("\t"));

    if (sample_line_split.size() != 3) {

        cerr << "\nERROR: Line \"" << sample_line << "\" in the samples file should contain three tab-seperated columns (<Sample ID>, <Gender> & <KMC Output Prefix>)\n" << endl;
        exit(1);
    }

    name = sample_line_split.at(0);

    if ((sample_line_split.at(1) == "F") or (sample_line_split.at(1) == "Female")) {

        gender = Utils::Gender::Female;

    } else {

        if ((sample_line_split.at(1) != "M") and (sample_line_split.at(1) != "Male")) {

            cerr << "\nERROR: Gender (column two) in line \"" << sample_line << "\" in the samples file should be either \"F\" (Female) or \"M\" (Male)\n" << endl;
            exit(1);
        }

        gender = Utils::Gender::Male;
    }

    file = sample_line_split.at(2);
}



