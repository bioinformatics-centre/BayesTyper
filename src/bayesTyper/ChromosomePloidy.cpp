
/*
ChromosomePloidy.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#include <iostream>
#include <fstream>

#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string.hpp"

#include "ChromosomePloidy.hpp"
#include "Utils.hpp"
#include "Sample.hpp"

ChromosomePloidy::ChromosomePloidy(const string & chrom_ploidy_filename, const Chromosomes & chromosomes, const vector<Sample> & samples) {

	if (chrom_ploidy_filename.empty()) {

		auto chromosomes_it = chromosomes.cbegin();

		while (chromosomes_it != chromosomes.cend()) {

			if (!chromosomes.isDecoy(chromosomes_it->first)) {

				auto gender_ploidy_it = gender_ploidy.emplace(chromosomes_it->first, vector<Utils::Ploidy>(2, Utils::Ploidy::Diploid));
				assert(gender_ploidy_it.second);

				auto sample_ploidy_it = sample_ploidy.emplace(chromosomes_it->first, vector<Utils::Ploidy>(samples.size(), Utils::Ploidy::Diploid));
				assert(sample_ploidy_it.second);

				auto chrom_name = chromosomes_it->first;		
		        transform(chrom_name.begin(), chrom_name.end(), chrom_name.begin(), ::tolower);

		        if ((chrom_name == "x") or (chrom_name == "chrx")) {

		        	gender_ploidy_it.first->second.at(static_cast<uchar>(Utils::Gender::Female)) = Utils::Ploidy::Diploid;
		        	gender_ploidy_it.first->second.at(static_cast<uchar>(Utils::Gender::Male)) = Utils::Ploidy::Haploid;

					for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

						if (samples.at(sample_idx).gender == Utils::Gender::Male) {

							sample_ploidy_it.first->second.at(sample_idx) = Utils::Ploidy::Haploid;					
						}
					} 

		        } else if ((chrom_name == "y") or (chrom_name == "chry")) {

		        	gender_ploidy_it.first->second.at(static_cast<uchar>(Utils::Gender::Female)) = Utils::Ploidy::Null;
		        	gender_ploidy_it.first->second.at(static_cast<uchar>(Utils::Gender::Male)) = Utils::Ploidy::Haploid;

					for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

						if (samples.at(sample_idx).gender == Utils::Gender::Female) {

							sample_ploidy_it.first->second.at(sample_idx) = Utils::Ploidy::Null;
						
						} else {

							sample_ploidy_it.first->second.at(sample_idx) = Utils::Ploidy::Haploid;
						}
					} 
		        }
	    	}

	        chromosomes_it++;
		}

	} else {

		ifstream chrom_ploidy_infile(chrom_ploidy_filename.c_str());

        if (!chrom_ploidy_infile.is_open()) {

            cerr << "\nERROR: Unable to open file " << chrom_ploidy_filename << "\n" << endl;
            exit(1);
        }

		unordered_map<string, pair<int, int> > chrom_ploidies;

		for (string chrom_ploidy_line; getline(chrom_ploidy_infile, chrom_ploidy_line);) {

			vector<string> chrom_ploidy_line_split;
		    boost::split(chrom_ploidy_line_split, chrom_ploidy_line, boost::is_any_of("\t"));

		    if (chrom_ploidy_line_split.size() != 3) {

		        cerr << "\nERROR: Line \"" << chrom_ploidy_line << "\" in the chromosome ploidy file should contain three tab-seperated columns (<Chromosome name>, <Female Ploidy> & <Male Ploidy>)\n" << endl;
		        exit(1);
		    }

			const int female_ploidy = stoi(chrom_ploidy_line_split.at(1));
			const int male_ploidy = stoi(chrom_ploidy_line_split.at(2));

		    if (!chrom_ploidies.emplace(chrom_ploidy_line_split.front(), make_pair(female_ploidy, male_ploidy)).second) {

		        cerr << "\nERROR: Chromosome (column one) in line \"" << chrom_ploidy_line << "\" appear multiple times in the chromosome ploidy file\n" << endl;
		        exit(1);
		    }

		    if ((female_ploidy < 0) or (female_ploidy > 2)) {

		        cerr << "\nERROR: Female ploidy (column two) in line \"" << chrom_ploidy_line << "\" in the chromosome ploidy file should be between zero and two; only ploidy levels up to diploid are currently supported\n" << endl;
		        exit(1);
		    }

		    if ((male_ploidy < 0) or (male_ploidy > 2)) {

		        cerr << "\nERROR: Male ploidy (column three) in line \"" << chrom_ploidy_line << "\" in the chromosome ploidy file should be between zero and two; only ploidy levels up to diploid are currently supported\n" << endl;
		        exit(1);
		    }
		}

	    chrom_ploidy_infile.close();

		auto chromosomes_it = chromosomes.cbegin();

		while (chromosomes_it != chromosomes.cend()) {

			if (!chromosomes.isDecoy(chromosomes_it->first)) {

				auto chrom_ploidies_it = chrom_ploidies.find(chromosomes_it->first);

			    if (chrom_ploidies_it == chrom_ploidies.end()) {

			        cerr << "\nERROR: Chromosome \"" << chromosomes_it->first << "\" in reference genome does not appear in the chromosome ploidy file\n" << endl;
			        exit(1);
			    }

				auto gender_ploidy_it = gender_ploidy.emplace(chrom_ploidies_it->first, vector<Utils::Ploidy>(2, Utils::Ploidy::Diploid));
				assert(gender_ploidy_it.second);

		        gender_ploidy_it.first->second.at(static_cast<uchar>(Utils::Gender::Female)) = static_cast<Utils::Ploidy>(chrom_ploidies_it->second.first);
		        gender_ploidy_it.first->second.at(static_cast<uchar>(Utils::Gender::Male)) = static_cast<Utils::Ploidy>(chrom_ploidies_it->second.second);

				auto sample_ploidy_it = sample_ploidy.emplace(chrom_ploidies_it->first, vector<Utils::Ploidy>(samples.size(), Utils::Ploidy::Diploid));
				assert(sample_ploidy_it.second);

				for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

					if (samples.at(sample_idx).gender == Utils::Gender::Female) {

						sample_ploidy_it.first->second.at(sample_idx) = static_cast<Utils::Ploidy>(chrom_ploidies_it->second.first);
					
					} else {

						sample_ploidy_it.first->second.at(sample_idx) = static_cast<Utils::Ploidy>(chrom_ploidies_it->second.second);
					}
				} 	
			}		

	        chromosomes_it++;
		}
	}
}

const vector<Utils::Ploidy> & ChromosomePloidy::getGenderPloidy(const string & chrom_name) const {

	auto gender_ploidy_it = gender_ploidy.find(chrom_name);
	assert(gender_ploidy_it != gender_ploidy.end());

	return gender_ploidy_it->second;
}

const vector<Utils::Ploidy> & ChromosomePloidy::getSamplePloidy(const string & chrom_name) const {

	auto sample_ploidy_it = sample_ploidy.find(chrom_name);
	assert(sample_ploidy_it != sample_ploidy.end());

	return sample_ploidy_it->second;
}
