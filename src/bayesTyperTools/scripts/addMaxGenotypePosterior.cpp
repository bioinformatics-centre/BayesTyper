
/*
addMaxGenotypePosterior.cpp - This file is part of BayesTyper (v0.9)


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
#include <numeric>

#include "VcfFile.hpp"
#include "Variant.hpp"
#include "JoiningString.hpp"
#include "Stats.hpp"
#include "Auxiliaries.hpp"


int main(int argc, char const *argv[]) {

	if (argc != 3) {

		std::cout << "USAGE: addMaxGenotypePosterior <input> <output_prefix>" << std::endl;
		return 1;
	}

    cout << "\n[" << Utils::getLocalTime() << "] Running BayesTyperTools (" << BT_VERSION << ") addMaxGenotypePosterior script ...\n" << endl;

	GenotypedVcfFileReader vcf_reader(string(argv[1]), true);

	auto output_meta_data = vcf_reader.metaData();
	output_meta_data.formatDescriptors().emplace("MaxGPP", Attribute::DetailedDescriptor("MaxGPP", Attribute::Number::One, Attribute::Type::Float, "Maximum genotype posterior probability"));

	VcfFileWriter vcf_writer(string(argv[2]) + ".vcf", output_meta_data, true);

	const vector<string> sample_ids = vcf_reader.metaData().sampleIds();
	vector<vector<uint> > max_allele_end_pos(sample_ids.size(), vector<uint>(2, 0));

	Variant * cur_var;

	ulong num_variants = 0;
	ulong num_gpp_genotypes = 0;
	
	while (vcf_reader.getNextVariant(&cur_var)) {

		num_variants++;

        for (ushort sample_idx = 0; sample_idx < sample_ids.size(); sample_idx++) {

            Sample * cur_sample = &(cur_var->getSample(sample_ids.at(sample_idx)));
			assert(cur_sample->ploidy() != Sample::Ploidy::Polyploid);
            
			auto max_gpp = Auxiliaries::getMaxGenotypePosterior(*cur_sample);

   			if (max_gpp.second) {

	        	num_gpp_genotypes++;
				cur_sample->info().setValue<float>("MaxGPP", max_gpp.first);
   			}
		}

        vcf_writer.write(cur_var);
		delete cur_var;

		if ((num_variants % 100000) == 0) {

			std::cout << "[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants" << endl;
		}
	}

	cout << "\n[" << Utils::getLocalTime() << "] Added maximum genotype posterior probability (MaxGPP) to " << num_gpp_genotypes << " genotypes." << endl;
	cout << endl;

	return 0;
}
