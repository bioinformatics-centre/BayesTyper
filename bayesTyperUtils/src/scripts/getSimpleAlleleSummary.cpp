
/*
getSimpleAlleleSummary.cpp - This file is part of BayesTyper (v0.9)


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

		std::cout << "USAGE: getSimpleAlleleSummary <vcf> <out_prefix>" << std::endl;
		return 1;
	}

	cout << "\n[" << Utils::getLocalTime() << "] Running BayesTyperUtils (" << BTU_VERSION << ") getSimpleAlleleSummary script ...\n" << endl;

	VcfFileReader vcf_reader(argv[1], false);
	Auxiliaries::removeNonRelevantInfoDescriptors(&(vcf_reader.metaData()), {"AC", "MLEA"});

	unordered_map<string, uint> allele_summary_stats;

	Variant * cur_var;
	uint num_variants = 0;

	while (vcf_reader.getNextVariant(&cur_var)) {

		num_variants++;

		for (uint alt_idx = 0; alt_idx < cur_var->numAlts(); alt_idx++) {

			auto allele_attributes = Auxiliaries::alleleAttributes(cur_var->alt(alt_idx), cur_var->ref());
			
			JoiningString allele_summary_line_elements('\t');
			allele_summary_line_elements.join(cur_var->chrom());

			if (cur_var->filters().empty()) {

				allele_summary_line_elements.join(".");

			} else {

				JoiningString filter_elements(';');
				filter_elements.join(cur_var->filters());

				allele_summary_line_elements.join(filter_elements.str());
			}

			allele_summary_line_elements.join(allele_attributes.typeStr());
			allele_summary_line_elements.join(to_string(allele_attributes.length));
			allele_summary_line_elements.join(to_string(allele_attributes.num_ambiguous));
			allele_summary_line_elements.join(to_string(allele_attributes.sv_length));

			auto ac_value = cur_var->alt(alt_idx).info().getValue<int>("AC");		
			auto mleac_value = cur_var->alt(alt_idx).info().getValue<int>("MLEAC");		

			if (ac_value.second) {

				allele_summary_line_elements.join(to_string(ac_value.first));
			
			} else if (mleac_value.second) {

				allele_summary_line_elements.join(to_string(mleac_value.first));
			
			} else {

				allele_summary_line_elements.join("-1");
			}

			auto allele_summary_stats_emplace = allele_summary_stats.emplace(allele_summary_line_elements.str(), 0);
			allele_summary_stats_emplace.first->second++;
		}
				
		if ((num_variants % 100000) == 0) {

			std::cout << "[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants" << endl;
		}

		delete cur_var;
	}

	ofstream output_writer(string(argv[2]) + ".txt");
	output_writer << "Count\tChromosome\tFilter\tAlleleType\tAlleleLength\tNumAmbiguous\tAlleleSVLength\tAlleleCount\n";

	for (auto &line: allele_summary_stats) {

		output_writer << line.second << "\t" << line.first << endl;
	}

	output_writer.close();

	cout << "\n[" << Utils::getLocalTime() << "] Wrote summary statistics\n" << endl;

	return 0;
}
