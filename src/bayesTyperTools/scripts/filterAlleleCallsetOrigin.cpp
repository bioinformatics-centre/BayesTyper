
/*
filterAlleleCallsetOrigin.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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

#include "VcfFile.hpp"
#include "Variant.hpp"
#include "JoiningString.hpp"
#include "Auxiliaries.hpp"


int main(int argc, char const *argv[]) {

	if ((argc != 4) and (argc != 5)) {

		std::cout << "USAGE: filterAlleleCallsetOrigin <variant_file> <output_prefix> <ACO_tag>,<ACO_tag>,... (<use_reference_allele>)" << std::endl;
		return 1;
	}

    cout << "\n[" << Utils::getLocalTime() << "] Running BayesTyperTools (" << BT_VERSION << ") filterAlleleCallsetOrigin script ...\n" << endl;
	
	GenotypedVcfFileReader vcf_reader(string(argv[1]), true);

	VcfFileWriter output_vcf(string(argv[2]) + ".vcf.gz", vcf_reader.metaData(), true);
	Variant * cur_var;

	auto aco_filter_tags = Utils::splitString(argv[3], ',');
	bool filter_missing = find(aco_filter_tags.begin(), aco_filter_tags.end(), ".") != aco_filter_tags.end();
	bool use_reference_allele = (argc == 5) ? stoi(argv[4]) : false;

	uint num_variants = 0;
	uint num_alt_alleles = 0;
	uint num_rm_variants = 0;
	uint num_rm_alt_alleles = 0;

	while (vcf_reader.getNextVariant(&cur_var)) {

		num_variants++;

		vector<uint> rm_alt_indices;
		rm_alt_indices.reserve(cur_var->numAlts());

		for (uint alt_idx = 0; alt_idx < cur_var->numAlts(); alt_idx++) {

			num_alt_alleles++;

            auto aco_value = cur_var->alt(alt_idx).info().getValue<string>("ACO"); 

            if (aco_value.second) {

				auto aco_value_split = Utils::splitString(aco_value.first, ':');
				uint num_filter_tags = 0;

				for (auto aco_filter_tag: aco_filter_tags) {

					if (find(aco_value_split.begin(), aco_value_split.end(), aco_filter_tag) != aco_value_split.end()) {
						
						num_filter_tags++;
					}
				}

				if (num_filter_tags == aco_value_split.size()) {

					rm_alt_indices.push_back(alt_idx);
				}

            } else if (filter_missing) {

            	rm_alt_indices.push_back(alt_idx);
            }
		}	

		num_rm_alt_alleles += rm_alt_indices.size();
			
		assert(rm_alt_indices.size() <= cur_var->numAlts());

		if (rm_alt_indices.empty()) {

			output_vcf.write(cur_var);

		} else {

			cur_var->removeAlts(rm_alt_indices, use_reference_allele);

			if (cur_var->numAlts() > 0) {

				if ((cur_var->numAlts() == 1) and (cur_var->alt(0).isMissing())) {

					num_rm_variants++;

				} else {

	            	Auxiliaries::updateAlleleStatsAndCallProb(cur_var);
					Auxiliaries::rightTrimVariant(cur_var);

					output_vcf.write(cur_var);
				}

			} else {

				num_rm_variants++;
			}
		}

		if ((num_variants % 100000) == 0) {

			std::cout << "[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants" << endl;
		}

		delete cur_var;
	}

	cout << "\n[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants and " << num_alt_alleles << " alternative alleles:" <<  endl;
    cout << "\n\t- Number of filtered variants: " << num_rm_variants << endl;
    cout << "\t- Number of filtered alternative alleles: " << num_rm_alt_alleles << endl;

	cout << endl;

	return 0;
}
