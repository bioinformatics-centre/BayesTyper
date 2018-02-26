
/*
getGenomicIntervals.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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

	if (argc != 4) {

		std::cout << "USAGE: getGenomicIntervals <input> <output_prefix> <min_var_per_interval>" << std::endl;
		return 1;
	}

    cout << "\n[" << Utils::getLocalTime() << "] Running BayesTyperTools (" << BT_VERSION << ") getGenomicIntervals script ...\n" << endl;
	
	VcfFileReader vcf_reader(string(argv[1]), true);

	vcf_reader.metaData().infoDescriptors().clear();
	vcf_reader.metaData().formatDescriptors().clear();
	
	ofstream interval_writer(string(argv[2]) + ".bed");
	assert(interval_writer.is_open());

	Variant * cur_var;
	uint num_variants = 0;

	uint min_var_per_interval = stoi(argv[3]);
	assert(min_var_per_interval > 0);

	uint num_intervals = 0;

	string cur_chrom = "";
	uint cur_chrom_len = 0;

	uint interval_num_variants = 0;
	uint interval_start_pos = 1;
	uint interval_end_pos = interval_start_pos;

	while (vcf_reader.getNextVariant(&cur_var)) {

		num_variants++;

		if (cur_var->chrom() != cur_chrom) {

			if (!(cur_chrom.empty())) {

				assert(cur_chrom_len > 0);
				assert(interval_start_pos <= cur_chrom_len);

				interval_writer << cur_chrom << "\t" << interval_start_pos - 1 << "\t" << cur_chrom_len - 1 << "\n";
				num_intervals++;
				
				interval_num_variants = 0;
				interval_start_pos = 1;
				interval_end_pos = interval_start_pos;
			}

			cur_chrom = cur_var->chrom();
			cur_chrom_len = vcf_reader.metaData().getContig(cur_var->chrom()).length();
		}

		if ((min_var_per_interval <= interval_num_variants) and (interval_end_pos < cur_var->pos())) {

			assert(interval_start_pos < cur_var->pos());

			interval_writer << cur_chrom << "\t" << interval_start_pos - 1 << "\t" << cur_var->pos() - 2 << "\n";
			num_intervals++;

			interval_num_variants = 0;
			interval_start_pos = cur_var->pos();	
			interval_end_pos = interval_start_pos;
		}

		interval_num_variants++;

		for (uint alt_allele_idx = 0; alt_allele_idx < cur_var->numAlts(); alt_allele_idx++) {

			assert(!(cur_var->alt(alt_allele_idx).isID()));

            if (cur_var->alt(alt_allele_idx).isMissing()) {

                continue;
            }

            Allele cur_ref_allele = cur_var->ref();
            Allele cur_alt_allele = cur_var->alt(alt_allele_idx);

            Auxiliaries::rightTrimAllelePair(&cur_ref_allele, &cur_alt_allele);

            assert(!(cur_ref_allele.seq().empty()));
            assert(!(cur_alt_allele.seq().empty()));

            interval_end_pos = max(interval_end_pos, static_cast<uint>(cur_var->pos() + cur_ref_allele.seq().size() - 1));
		}	

		if ((num_variants % 100000) == 0) {

			std::cout << "[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants" << endl;
		}

		delete cur_var;
	}

	assert(!(cur_chrom.empty()));
	assert(cur_chrom_len > 0);
	assert(interval_start_pos <= cur_chrom_len);

	interval_writer << cur_chrom << "\t" << interval_start_pos - 1 << "\t" << cur_chrom_len - 1 << "\n";
	num_intervals++;

	interval_writer.close();	

	cout << "\n[" << Utils::getLocalTime() << "] Wrote " << num_intervals << " genomic intervals" << endl;
	cout << endl;

	return 0;
}
