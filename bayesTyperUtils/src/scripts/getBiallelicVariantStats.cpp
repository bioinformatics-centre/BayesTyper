
/*
getBiallelicVariantStats.cpp - This file is part of BayesTyper (v0.9)


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
#include <unordered_map>

#include "vcf++/VcfFile.hpp"
#include "vcf++/Variant.hpp"
#include "vcf++/JoiningString.hpp"
#include "vcf++/Auxiliaries.hpp"

int main(int argc, char const *argv[]) {

	if (argc != 3) {

		std::cout << "USAGE: getBiallelicVariantStats <vcf> <out_prefix>" << std::endl;
		return 1;
	}

	cout << "\n[" << Utils::getLocalTime() << "] Running BayesTyperUtils (" << BTU_VERSION << ") biallelicVariantStats script ...\n" << endl;

	string vcf_filename(argv[1]);
	VcfFileReader vcf_reader(vcf_filename, true);
	Auxiliaries::removeNonRelevantInfoDescriptors(&(vcf_reader.metaData()), {"AC", "AN", "AT"});

	ofstream output_writer(string(argv[2]) + ".txt");
	output_writer << "Chr\tPos\tType\tMaxLength\tMaxTrimLength\tAlleleCount\tAlleleNumber\n";

	Variant * cur_var;
	uint num_variants = 0;

	while (vcf_reader.getNextVariant(&cur_var)) {

		num_variants++;

		if (cur_var->numAlls() == 2) {

			assert(!(cur_var->ref().isMissing()));
			assert(!(cur_var->allele(1).isMissing()));

			uint max_length = max(cur_var->ref().seq().size(), cur_var->allele(1).seq().size());

			auto allele_attributes = Auxiliaries::alleleAttributes(cur_var->allele(1), cur_var->ref());
			auto at_value = cur_var->allele(1).info().getValue<string>("AT");
			
			assert(at_value.second);
			assert(allele_attributes.typeStr() == at_value.first);

			uint num_trimmed = Auxiliaries::rightTrimAllelePair(&(cur_var->ref()), &(cur_var->allele(1)));
			uint max_trim_length = max(cur_var->ref().seq().size(), cur_var->allele(1).seq().size());

			assert(max_length == (max_trim_length + num_trimmed));

			auto ac_value = cur_var->allele(1).info().getValue<int>("AC");		
			assert(ac_value.second);			

			auto an_value = cur_var->info().getValue<int>("AN");		
			assert(an_value.second);	
			
			JoiningString stats_line_elements('\t');
			stats_line_elements.join(cur_var->chrom());
			stats_line_elements.join(to_string(cur_var->pos()));
			stats_line_elements.join(at_value.first);
			stats_line_elements.join(to_string(max_length));
			stats_line_elements.join(to_string(max_trim_length));
			stats_line_elements.join(to_string(ac_value.first));
			stats_line_elements.join(to_string(an_value.first));

			output_writer << stats_line_elements.str() << endl;
		}
				
		if ((num_variants % 100000) == 0) {

			std::cout << "[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants" << endl;
		}

		delete cur_var;
	}

	cout << "\n[" << Utils::getLocalTime() << "] Wrote biallelic variant stats statistics" << endl;

	return 0;
}
