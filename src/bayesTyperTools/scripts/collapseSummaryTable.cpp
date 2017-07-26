
/*
collapseSummaryTable.cpp - This file is part of BayesTyper (v0.9)


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

		std::cout << "USAGE: collapseSummaryTable <input> <output_prefix> <num_count_columns (assume first columns)>" << std::endl;
		return 1;
	}

    cout << "\n[" << Utils::getLocalTime() << "] Running BayesTyperTools (" << BT_VERSION << ") collapseSummaryTable script ...\n" << endl;

    ifstream reader(argv[1]);

    auto num_count_columns = stoi(argv[3]);
    assert(num_count_columns > 0);

    string line;
    uint line_count = 0;

    string header;
    unordered_map<string, vector<ulong> > collapsed_lines;

    while (getline(reader, line)) {

        line_count++;

        if (line_count == 1) {

            assert(header.empty());
            header = line;

            assert(num_count_columns <= (count(line.begin(), line.end(), '\t')));
            continue;
        }

        auto line_split = Utils::splitString(line, '\t', num_count_columns + 1);
        assert(line_split.size() == static_cast<uint>(num_count_columns + 1));

        auto collapsed_lines_it = collapsed_lines.emplace(line_split.back(), vector<ulong>(num_count_columns, 0));

        for (uint count_idx = 0; count_idx < (line_split.size() - 1); count_idx++) {

            collapsed_lines_it.first->second.at(count_idx) += stol(line_split.at(count_idx));
        }

        if ((line_count % 1000000) == 0) {

            std::cout << "[" << Utils::getLocalTime() << "] Parsed " << line_count << " lines" << endl;
        }
    }

    reader.close();

    assert(!(header.empty()));

    cout << "\n[" << Utils::getLocalTime() << "] Writing collapsed file ..." << endl;

    ofstream writer(string(argv[2]) + ".txt");
    writer << header << "\n";

    for (auto & line: collapsed_lines) {

        for (auto & count: line.second) {

            writer << count << "\t";
        }

        writer << line.first << "\n";
    }

    writer.close();

	cout << "[" << Utils::getLocalTime() << "] Collapsed " << line_count << " lines down to " << collapsed_lines.size() + 1 << endl;
	cout << endl;

	return 0;
}
