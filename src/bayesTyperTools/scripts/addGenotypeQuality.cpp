
/*
addGenotypeQuality.cpp - This file is part of BayesTyper (v0.9)


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

		std::cout << "USAGE: addGenotypeQuality <input> <output_prefix>" << std::endl;
		return 1;
	}

    cout << "\n[" << Utils::getLocalTime() << "] Running BayesTyperTools (" << BT_VERSION << ") addGenotypeQuality script ...\n" << endl;

	GenotypedVcfFileReader vcf_reader(argv[1], true);

	auto output_meta_data = vcf_reader.metaData();
	output_meta_data.formatDescriptors().emplace("GQ", Attribute::DetailedDescriptor({make_pair("ID","GQ"), make_pair("Number","1"), make_pair("Type","Integer"), make_pair("Description","Genotype quality (estimated from GL)")}));

	VcfFileWriter vcf_writer(string(argv[2]) + ".vcf", output_meta_data, true);
	
	Variant * cur_var;
	int num_variants = 0;

	while (vcf_reader.getNextVariant(&cur_var)) {

		num_variants++;

    	for (auto & sample_id: vcf_reader.metaData().sampleIds()) {

        	Sample * cur_sample = &(cur_var->getSample(sample_id));

        	assert(cur_sample->ploidy() != Sample::Ploidy::Polyploid);

        	if (cur_sample->ploidy() != Sample::Ploidy::Zeroploid) {

        		if (cur_sample->genotypeInfo().front().getValue<float>("GL").second) {

		        	vector<float> phred_likelihoods;
		        	phred_likelihoods.reserve(((cur_var->numAlls() + 1) * cur_var->numAlls()) / 2);

		        	for (auto & attribute: cur_sample->genotypeInfo()) {

		        		auto gl_value = attribute.getValue<float>("GL");
		        		assert(gl_value.second);

		        		phred_likelihoods.push_back(-10 * gl_value.first);
		        	}

		        	assert(phred_likelihoods.size() >= 3);
		        	sort(phred_likelihoods.begin(), phred_likelihoods.end());
		        	
		        	if (!(Utils::floatCompare(phred_likelihoods.at(0), 0))) {

		        		assert((cur_sample->callStatus() == Sample::CallStatus::Missing) or (cur_sample->callStatus() == Sample::CallStatus::Partial));
		        	}

					assert(cur_sample->info().setValue("GQ", round(phred_likelihoods.at(1) - phred_likelihoods.at(0))));
				}
			}
        }

		vcf_writer.write(cur_var);
		delete cur_var;

		if ((num_variants % 100000) == 0) {

			std::cout << "[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants" << endl;
		}		
	}

	cout << "\n[" << Utils::getLocalTime() << "] Added genotype quality to " << vcf_reader.metaData().numSamples() << " samples across " << num_variants << " variants" << endl;
	cout << endl;

	return 0;
}
