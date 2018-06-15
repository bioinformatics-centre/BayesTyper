
/*
filterStructuralVariants.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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

	if (argc != 6) {

		std::cout << "USAGE: filterStructuralVariants <variant_file> <output_prefix> <max_allele_length (reference | alternative)> <min_sv_length (alternative_length - reference_length)> <max_sv_length (alternative_length - reference_length)>" << std::endl;
		return 1;
	}

    cout << "\n[" << Utils::getLocalTime() << "] Running BayesTyperTools (" << BT_VERSION << ") filterStructuralVariants script ...\n" << endl;
	
	VcfFileReader vcf_reader(string(argv[1]), true);

	VcfFileWriter output_vcf(string(argv[2]) + ".vcf.gz", vcf_reader.metaData(), true);
	Variant * cur_var;

	uint max_allele_length = stoi(argv[3]);
	int min_sv_length = stoi(argv[4]);
	int max_sv_length = stoi(argv[5]);

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

			auto allele_attributes = Auxiliaries::alleleAttributes(cur_var->alt(alt_idx), cur_var->ref());

			assert(static_cast<long>(allele_attributes.length) >= allele_attributes.sv_length);

			uint ref_length = allele_attributes.length - allele_attributes.sv_length;
			uint alt_length = allele_attributes.length;

			if ((max_allele_length < ref_length) or (max_allele_length < alt_length)) {

				rm_alt_indices.push_back(alt_idx);
			
			} else if ((allele_attributes.sv_length < min_sv_length) or (max_sv_length < allele_attributes.sv_length)) {

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
	        cur_var->setFilters({});

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
    cout << "\n\t- Number of filtered structural variants: " << num_rm_variants << endl;
    cout << "\t- Number of filtered alternative alleles: " << num_rm_alt_alleles << endl;

	cout << endl;

	return 0;
}
