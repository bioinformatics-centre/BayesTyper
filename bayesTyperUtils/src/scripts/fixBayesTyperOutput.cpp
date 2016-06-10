
/*
fixBayesTyperOutput.cpp - This file is part of BayesTyper (v0.9)


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

void convertAlleleDescriptorToDot(VcfMetaData * meta_data, const string & descriptor_name) {

	auto descriptor = meta_data->infoDescriptors().find(descriptor_name);
	
	if (descriptor != meta_data->infoDescriptors().end()) {

		vector<pair<string, string> > new_destriptor_elements({make_pair("ID",descriptor->second.id()), make_pair("Number","."), make_pair("Type","String"), make_pair("Description",descriptor->second.description())});
		descriptor->second = Attribute::DetailedDescriptor(new_destriptor_elements);
	}
}

void convertDotAttributeToAllele(Variant * cur_var, const string & descriptor_name) {

	auto descriptor_var = cur_var->info().getValue(descriptor_name);

	if (descriptor_var.second) {

		auto descriptor_var_split = Utils::splitString(descriptor_var.first.str(), ',');

		if (descriptor_var_split.size() != cur_var->numAlts()) {

			if (cur_var->alt(cur_var->numAlts() - 1).isMissing()) {

				descriptor_var_split.push_back(".");
			}
		}

		assert(descriptor_var_split.size() == cur_var->numAlts());

		for (uint i = 0; i < cur_var->numAlts(); i++) {

			assert(cur_var->alt(i).info().setValue<string>(descriptor_name, descriptor_var_split.at(i)));
		}
	}
}

int main(int argc, char const *argv[]) {

	if (argc != 3) {

		std::cout << "USAGE: fixBayesTyperOutput <bayesTyper_output> <output_prefix>" << std::endl;
		return 1;
	}

    cout << "\n[" << Utils::getLocalTime() << "] Running BayesTyperUtils (" << BTU_VERSION << ") fixBayesTyperOutput script ...\n" << endl;

	string vcf_filename(argv[1]);
	GenotypedVcfFileReader vcf_file(vcf_filename, true);

	VcfFileWriter output_vcf(string(argv[2]) + ".vcf", vcf_file.metaData(), true);
	Variant * current_variant;

	convertAlleleDescriptorToDot(&(vcf_file.metaData()), "ACO");
	convertAlleleDescriptorToDot(&(vcf_file.metaData()), "AsmVar_ASQR");

	int vars = 0;

	while (vcf_file.getNextVariant(&current_variant)) {

		convertDotAttributeToAllele(current_variant, "ACO");
		convertDotAttributeToAllele(current_variant, "AsmVar_ASQR");

		vars++;
		output_vcf.write(current_variant);

		if ((vars % 100000) == 0) {

			std::cout << "[" << Utils::getLocalTime() << "] Parsed " << vars << " variants" << endl;
		}

		delete current_variant;
	}

	cout << "\n[" << Utils::getLocalTime() << "] Fixed " << vars << " BayesTyper variants" << endl;
	cout << endl;

	return 0;
}
