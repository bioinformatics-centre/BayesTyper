
/*
collapseSummaryTable.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include "assert.h"

#include "Utils.hpp"
#include "JoiningString.hpp"

using namespace std;

int main(int argc, char const *argv[]) {

	if (argc != 4) {

		std::cout << "USAGE: collapseSummaryTable <variant_file> <output_prefix> <num_count_columns (assume first columns)>" << std::endl;
		return 1;
	}

    cout << "\n[" << Utils::getLocalTime() << "] Running BayesTyperTools (" << BT_VERSION << ") collapseSummaryTable script ...\n" << endl;

    ifstream summary_infile(argv[1]);

    if (!summary_infile.is_open()) {

        cerr << "\nERROR: Unable to open file " << argv[1] << "\n" << endl;
        exit(1);
    }

    auto num_count_columns = stoi(argv[3]);
    assert(num_count_columns > 0);

    uint summary_line_count = 0;

    string header;
    unordered_map<string, vector<ulong> > collapsed_lines;

    for (string summary_line; getline(summary_infile, summary_line);) {

        summary_line_count++;

        if (summary_line_count == 1) {

            assert(header.empty());
            header = summary_line;

            assert(num_count_columns <= (count(summary_line.begin(), summary_line.end(), '\t')));
            continue;
        }

        auto summary_line_split = Utils::splitString(summary_line, '\t', num_count_columns + 1);
        assert(summary_line_split.size() == static_cast<uint>(num_count_columns + 1));

        auto collapsed_lines_it = collapsed_lines.emplace(summary_line_split.back(), vector<ulong>(num_count_columns, 0));

        for (uint count_idx = 0; count_idx < (summary_line_split.size() - 1); count_idx++) {

            collapsed_lines_it.first->second.at(count_idx) += stol(summary_line_split.at(count_idx));
        }

        if ((summary_line_count % 1000000) == 0) {

            std::cout << "[" << Utils::getLocalTime() << "] Parsed " << summary_line_count << " lines" << endl;
        }
    }

    summary_infile.close();

    assert(!header.empty());

    cout << "\n[" << Utils::getLocalTime() << "] Writing collapsed file ..." << endl;

    ofstream summary_outfile(string(argv[2]) + ".txt");

    if (!summary_outfile.is_open()) {

        cerr << "\nERROR: Unable to write file " << string(argv[2]) + ".txt" << "\n" << endl;
        exit(1);
    }

    summary_outfile << header << "\n";

    for (auto & line: collapsed_lines) {

        for (auto & count: line.second) {

            summary_outfile << count << "\t";
        }

        summary_outfile << line.first << "\n";
    }

    summary_outfile.close();

	cout << "[" << Utils::getLocalTime() << "] Collapsed " << summary_line_count << " lines down to " << collapsed_lines.size() + 1 << endl;
	cout << endl;

	return 0;
}
