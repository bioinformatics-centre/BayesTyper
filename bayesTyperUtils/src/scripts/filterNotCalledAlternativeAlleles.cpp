
/*
filterNotCalledAlternativeAlleles.cpp - This file is part of BayesTyper (v0.9)


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

#include "vcf++/VcfFile.hpp"
#include "vcf++/Variant.hpp"
#include "vcf++/JoiningString.hpp"
#include "vcf++/Stats.hpp"
#include "vcf++/Auxiliaries.hpp"


int main(int argc, char const *argv[]) {

	if (argc != 3) {

		std::cout << "USAGE: filterNotCalledAlternativeAlleles <variants> <output_prefix>" << std::endl;
		return 1;
	}

    cout << "\n[" << Utils::getLocalTime() << "] Running BayesTyperUtils (" << BTU_VERSION << ") filterNotCalledAlternativeAlleles script ...\n" << endl;
	
	GenotypedVcfFileReader vcf_reader(string(argv[1]), true);
	Auxiliaries::removeNonRelevantFormatDescriptors(&(vcf_reader.metaData()), {"GT"});

	VcfFileWriter output_vcf(string(argv[2]) + ".vcf", vcf_reader.metaData(), true);
	Variant * cur_var;

	uint num_variants = 0;
	uint num_alt_alleles = 0;
	uint num_rm_variants = 0;
	uint num_rm_alt_alleles = 0;

	while (vcf_reader.getNextVariant(&cur_var)) {

		num_variants++;

		auto allele_stats = Stats::calcAlleleStats(cur_var);
		assert(allele_stats.first.allele_counts.size() > 1);
		
		assert(allele_stats.first.allele_counts.size() == cur_var->numAlls());
		assert(allele_stats.first.allele_counts.size() == (cur_var->numAlts() + 1));

		vector<uint> rm_alt_indices;
		rm_alt_indices.reserve(cur_var->numAlts());

		for (uint alt_idx = 0; alt_idx < cur_var->numAlts(); alt_idx++) {

			num_alt_alleles++;

			if (allele_stats.first.allele_counts.at(alt_idx + 1) == 0) {

				num_rm_alt_alleles++;
				rm_alt_indices.push_back(alt_idx);
			}
		}	

		num_rm_alt_alleles += rm_alt_indices.size();
			
		assert(rm_alt_indices.size() <= cur_var->numAlts());
		assert(rm_alt_indices.back() < cur_var->numAlts());

		cur_var->removeAlts(rm_alt_indices);

		if (cur_var->numAlts() > 0) {

	        cur_var->setIds({});
	        cur_var->setQual({0, false});
	        cur_var->setFilters({"PASS"});
		
			Auxiliaries::rightTrimVariant(cur_var);
			output_vcf.write(cur_var);

		} else {

			num_rm_variants++;
		}

		if ((num_variants % 100000) == 0) {

			std::cout << "[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants" << endl;
		}

		delete cur_var;
	}

	cout << "\n[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants and " << num_alt_alleles << " alternative alleles:" <<  endl;
    cout << "\n\t- Number of not called variants filtered: " << num_rm_variants << endl;
    cout << "\t- Number of not called alternative alleles filtered: " << num_rm_alt_alleles << endl;

	cout << endl;

	return 0;
}
