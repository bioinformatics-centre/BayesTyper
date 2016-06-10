
/*
addUniqueIdToVariants.cpp - This file is part of BayesTyper (v0.9)


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

#include "vcf++/VcfFile.hpp"
#include "vcf++/Variant.hpp"
#include "vcf++/JoiningString.hpp"



int main(int argc, char const *argv[]) {

	if (argc != 4) {

		std::cout << "USAGE: addUniqueIdToVariants <variants> <output_prefix> <id_prefix>" << std::endl;
		return 1;
	}

    cout << "\n[" << Utils::getLocalTime() << "] Running BayesTyperUtils (" << BTU_VERSION << ") addUniqueIdToVariants script ...\n" << endl;

	GenotypedVcfFileReader vcf_reader(string(argv[1]), true);
	VcfFileWriter output_vcf(string(argv[2]) + ".vcf", vcf_reader.metaData(), true);

	const string id_prefix = argv[3];

	Variant * cur_var;
	uint num_variants = 0;

	while (vcf_reader.getNextVariant(&cur_var)) {

		num_variants++;
		assert(num_variants <= 999999999);

		string cur_id = to_string(num_variants);
		cur_id.reserve(id_prefix.size() + 9);

		while (cur_id.size() < 9) {

			cur_id.insert(0, "0");
		}

		cur_id.insert(0, id_prefix);

    	cur_var->setIds({cur_id});
		output_vcf.write(cur_var);

		if ((num_variants % 100000) == 0) {

			std::cout << "[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants" << endl;
		}

		delete cur_var;
	}

	cout << "\n[" << Utils::getLocalTime() << "] Added ID's to " << num_variants << " variants." << endl;

	cout << endl;

	return 0;
}
