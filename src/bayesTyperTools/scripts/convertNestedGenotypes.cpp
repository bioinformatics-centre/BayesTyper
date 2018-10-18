
/*
convertNestedGenotypes.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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

		std::cout << "USAGE: convertNestedGenotypes <variant_file> <output_prefix>" << std::endl;
		return 1;
	}

    cout << "\n[" << Utils::getLocalTime() << "] Running BayesTyperTools (" << BT_VERSION << ") convertNestedGenotypes script ...\n" << endl;

	GenotypedVcfFileReader vcf_reader(string(argv[1]), true);
	VcfFileWriter vcf_writer(string(argv[2]) + ".vcf.gz", vcf_reader.metaData(), true);

	const vector<string> sample_ids = vcf_reader.metaData().sampleIds();
	vector<vector<uint> > max_allele_end_pos(sample_ids.size(), vector<uint>(2, 0));

	Variant * cur_var;
	string cur_chrom = "";

	ulong num_variants = 0;
	ulong num_genotypes = 0;
	ulong num_converted_genotypes = 0;
	
	while (vcf_reader.getNextVariant(&cur_var)) {

		if (cur_var->chrom() != cur_chrom) {

			max_allele_end_pos = vector<vector<uint> >(sample_ids.size(), vector<uint>(2, 0));
			cur_chrom = cur_var->chrom();
		}

		num_variants++;
		vector<uint> allele_ref_lengths(cur_var->numAlls(), 1);

		for (uint alt_allele_idx = 0; alt_allele_idx < cur_var->numAlts(); alt_allele_idx++) {

        	if (!cur_var->alt(alt_allele_idx).isMissing()) {

	        	Allele ref_allele = cur_var->ref();
	        	Allele alt_allele = cur_var->alt(alt_allele_idx);

	            Auxiliaries::rightTrimAllelePair(&ref_allele, &alt_allele);

	            assert(!ref_allele.seq().empty());
	            assert(!alt_allele.seq().empty());

	            allele_ref_lengths.at(alt_allele_idx + 1) = ref_allele.seq().size();
        	}
        }

        for (ushort sample_idx = 0; sample_idx < sample_ids.size(); sample_idx++) {

        	num_genotypes++;
            Sample * cur_sample = &(cur_var->getSample(sample_ids.at(sample_idx)));

			assert(cur_sample->ploidy() != Sample::Ploidy::Polyploid);
            
            const vector<ushort> genotype_estimate = cur_sample->genotypeEstimate();
            assert(genotype_estimate.size() <= 2);

            vector<ushort> new_genotype_estimate = genotype_estimate;

            for (uint genotyped_allele_idx = 0; genotyped_allele_idx < genotype_estimate.size(); genotyped_allele_idx++) {

            	if (cur_sample->isPhased()) {

					if (cur_var->pos() <= max_allele_end_pos.at(sample_idx).at(genotyped_allele_idx)) {

                		new_genotype_estimate.at(genotyped_allele_idx) = 0;
					}

            	} else if (cur_var->allele(genotype_estimate.at(genotyped_allele_idx)).isMissing()) {

            		new_genotype_estimate.at(genotyped_allele_idx) = 0;
            	}

				if (max_allele_end_pos.at(sample_idx).at(genotyped_allele_idx) < cur_var->pos()) {

	            	max_allele_end_pos.at(sample_idx).at(genotyped_allele_idx) = max(max_allele_end_pos.at(sample_idx).at(genotyped_allele_idx), static_cast<uint>(cur_var->pos() + allele_ref_lengths.at(genotype_estimate.at(genotyped_allele_idx)) - 1));
            	}
            }

            if (genotype_estimate != new_genotype_estimate) {

            	num_converted_genotypes++;
            	cur_sample->newGenotypeEstimate(new_genotype_estimate);
            }     	
		}

        vcf_writer.write(cur_var);
		delete cur_var;

		if ((num_variants % 100000) == 0) {

			std::cout << "[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants" << endl;
		}
	}

	cout << "\n[" << Utils::getLocalTime() << "] Out of " << num_genotypes << " genotypes in " << num_variants << " variants " << num_converted_genotypes << " nested genotypes were converted to the reference allele." << endl;
	cout << endl;

	return 0;
}
