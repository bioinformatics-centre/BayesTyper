
/*
getKmerStats.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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
#include <assert.h>
#include <limits>
#include <iomanip>

#include "kmc_api/kmc_file.h"

#include "Utils.hpp"
#include "JoiningString.hpp"

typedef unsigned char uchar;


static const uchar uchar_overflow = numeric_limits<uchar>::max();
static const uint kmer_size = BT_KMER_SIZE;


int main(int argc, char const *argv[]) {

	if (argc != 3) {

		std::cout << "USAGE: getKmerStats <kmc_table_prefix> <output_prefix>" << std::endl;
		return 1;
	}

    cout << "\n[" << Utils::getLocalTime() << "] Running BayesTyperTools (" << BT_VERSION << ") getKmerStats script ...\n" << endl;


	CKMCFile kmer_database;

	if (!kmer_database.OpenForListing(argv[1])) {

		cerr << "\nERROR: Unable to open KMC table " << argv[1] << "\n" << endl;
		exit(1);
	}

	uint32 db_kmer_size, mode, counter_size, lut_prefix_length, signature_len, min_count;
	uint64 max_count, total_kmers;
	
	assert(kmer_database.Info(db_kmer_size, mode, counter_size, lut_prefix_length, signature_len, min_count, max_count, total_kmers));
	assert(!mode);

	assert(db_kmer_size == kmer_size);

	cout << "[" << Utils::getLocalTime() << "] Parsing kmer table containing " << total_kmers << " unique kmers with a length of " << db_kmer_size << " nts ...\n" << endl;


	ulong num_kmers = 0;

	unordered_map<string, ulong> kmer_stats;
	
	CKmerAPI kmer_object(db_kmer_size);
	uint32 kmer_count;

	while (kmer_database.ReadNextKmer(kmer_object, kmer_count)) {

		num_kmers++;
		vector<ulong> nucleotide_counts(4, 0);

		for (auto & nucleotide: kmer_object.to_string()) {

			if (nucleotide == 'A') {

				nucleotide_counts.at(0)++;

			} else if (nucleotide == 'C') {

				nucleotide_counts.at(1)++;

			} else if (nucleotide == 'G') {

				nucleotide_counts.at(2)++;

			} else {

				assert(nucleotide == 'T');
				nucleotide_counts.at(3)++;
			}
		}

		assert(kmer_count <= uchar_overflow);

		JoiningString kmer_stats_line_elements('\t');
		kmer_stats_line_elements.join(to_string(kmer_count));
		kmer_stats_line_elements.join(to_string(nucleotide_counts.at(0)));
		kmer_stats_line_elements.join(to_string(nucleotide_counts.at(1)));
		kmer_stats_line_elements.join(to_string(nucleotide_counts.at(2)));
		kmer_stats_line_elements.join(to_string(nucleotide_counts.at(3)));

		auto kmer_stats_emplace = kmer_stats.emplace(kmer_stats_line_elements.str(), 0);
		kmer_stats_emplace.first->second++;

		if ((num_kmers % 10000000) == 0) {

			std::cout << "[" << Utils::getLocalTime() << "] Parsed " << num_kmers << " kmers" << endl;
		}		
	}

    ofstream stats_outfile(string(argv[2]) + "_kmer_stats.txt");

    if (!stats_outfile.is_open()) {

        cerr << "\nERROR: Unable to write file " << string(argv[2]) + "_kmer_stats.txt" << "\n" << endl;
        exit(1);
    }
    
	stats_outfile << "NumberOfKmers\tKmerCount\tAdenineCount\tCytosineCount\tGuanineCount\tThymineCount\n";

	for (auto & kmer_stat: kmer_stats) {

		stats_outfile << kmer_stat.second << "\t" << kmer_stat.first << "\n";
	}

	stats_outfile.close();

	cout << "\n[" << Utils::getLocalTime() << "] Wrote statistics for " << num_kmers << " kmers" << endl;
	cout << endl;

	return 0;
}
