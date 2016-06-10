
/*
splitMultiAllelicVariants.cpp - This file is part of BayesTyper (v0.9)


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
#include "vcf++/Stats.hpp"
#include "vcf++/Auxiliaries.hpp"


void updateGenotypes(Variant * cur_var, const vector<string> & complete_sample_ids) {

	for (auto & sample_id: complete_sample_ids) {

		Sample & cur_sample = cur_var->getSample(sample_id);
		vector<ushort> new_genotype_estimate = cur_sample.genotypeEstimate();

		if (cur_sample.ploidy() == Sample::Ploidy::Diploid) {

			assert(new_genotype_estimate.size() <= 2);

			if (new_genotype_estimate.size() < 2) {

				while (new_genotype_estimate.size() < 2) {

					new_genotype_estimate.insert(new_genotype_estimate.begin(), 0);
				}

				cur_sample.newGenotypeEstimate(new_genotype_estimate);
			}
 		
 		} else if (cur_sample.ploidy() == Sample::Ploidy::Haploid) {

			assert(new_genotype_estimate.size() <= 1);

			if (new_genotype_estimate.size() < 1) {

				new_genotype_estimate.insert(new_genotype_estimate.begin(), 0);
				cur_sample.newGenotypeEstimate(new_genotype_estimate);
			}
 		}
	}
}


void writeVariant(VcfFileWriter * output_vcf, Variant * cur_var) {

    cur_var->setIds({});
    cur_var->setQual({0, false});
    cur_var->setFilters({"PASS"});

	Auxiliaries::rightTrimVariant(cur_var);
    auto allele_stats = Stats::calcAlleleStats(cur_var);

    assert(allele_stats.first.allele_counts.size() == 2);
    assert(allele_stats.first.allele_freqs.size() == 2);

    assert(!cur_var->info().setValue<int>("AN", allele_stats.first.allele_count_sum));
   	assert(!cur_var->alt(0).info().setValue<int>("AC", allele_stats.first.allele_counts.at(1)));
    assert(!cur_var->alt(0).info().setValue<float>("AF", allele_stats.first.allele_freqs.at(1)));

	output_vcf->write(cur_var);
}


int main(int argc, char const *argv[]) {

	if (argc != 5) {

		std::cout << "USAGE: splitMultiAllelicVariants <variants> <output_prefix> <min_allele_posterior> <min_genotype_posterior>" << std::endl;
		return 1;
	}

    cout << "\n[" << Utils::getLocalTime() << "] Running BayesTyperUtils (" << BTU_VERSION << ") splitMultiAllelicVariants script ...\n" << endl;

	GenotypedVcfFileReader vcf_reader(string(argv[1]), true);
 	Auxiliaries::removeNonRelevantFormatDescriptors(&(vcf_reader.metaData()), {"GT", "GPP"});

	auto output_meta_data = vcf_reader.metaData();

	output_meta_data.filterDescriptors().clear();
	output_meta_data.infoDescriptors().erase("ACP");
	output_meta_data.infoDescriptors().erase("AsmVar_ASQR");
	output_meta_data.infoDescriptors().erase("DNE");
	output_meta_data.formatDescriptors().erase("GPP");

	VcfFileWriter output_vcf(string(argv[2]) + ".vcf", output_meta_data, true);

	float min_allele_posterior = stof(argv[3]);
	float min_genotype_posterior = stof(argv[4]);

	Variant * cur_var;
	auto sample_ids = output_meta_data.sampleIds();

	uint num_variants = 0;
	uint num_alt_alleles = 0;
	uint num_missing_alleles = 0;
	uint num_called_split_variants = 0;

	uint num_filtered_alleles = 0;
	uint num_filtered_samples = 0;

	while (vcf_reader.getNextVariant(&cur_var)) {

		num_variants++;
		num_alt_alleles += cur_var->numAlts();

		vector<string> complete_sample_ids;
		complete_sample_ids.reserve(sample_ids.size());

    	for (auto & sample_id: sample_ids) {

			Sample & cur_sample = cur_var->getSample(sample_id);

            if (cur_sample.callStatus() != Sample::CallStatus::Missing) {

            	assert(cur_sample.callStatus() == Sample::CallStatus::Complete);

            	auto genotype_gpp = Auxiliaries::getGenotypePosterior(cur_sample);
            	assert(genotype_gpp == Auxiliaries::getMaxGenotypePosterior(cur_sample));

            	if (genotype_gpp.second) {

		            if (genotype_gpp.first < min_genotype_posterior) {

		                cur_sample.clearGenotypeEstimate();
		            	assert(cur_sample.callStatus() == Sample::CallStatus::Missing);

		                num_filtered_samples++;
		            
		            } else {

		            	complete_sample_ids.push_back(sample_id);
		            }

            	} else {

            		assert(cur_sample.ploidy() == Sample::Ploidy::Zeroploid);
            	}
            } 
        }
    	
		while (cur_var->numAlts() > 1) {

			assert(!(cur_var->alt(0).isMissing()));

			auto acp_value = cur_var->alt(0).info().getValue<float>("ACP");
			assert(acp_value.second);

			if (acp_value.first < min_allele_posterior) {

				num_filtered_alleles++;

			} else {

				Variant * new_var = new Variant(*cur_var);

				vector<uint> alts_remove(new_var->numAlts() - 1);
				iota(alts_remove.begin(), alts_remove.end(), 1);

				new_var->removeAlts(alts_remove);
	
				assert(new_var->numAlls() == 2);
				assert(new_var->numAlts() == 1);

				num_called_split_variants++;

				updateGenotypes(new_var, complete_sample_ids);
				writeVariant(&output_vcf, new_var);
			}

			cur_var->removeAlts({0});
		}

		assert(cur_var->numAlls() == 2);
		assert(cur_var->numAlts() == 1);

		auto acp_value = cur_var->alt(0).info().getValue<float>("ACP");
		assert(acp_value.second);

		if (cur_var->alt(0).isMissing()) {

			num_missing_alleles++;
			
		} else if (acp_value.first < min_allele_posterior) {

			num_filtered_alleles++;

		} else {

			num_called_split_variants++;

			updateGenotypes(cur_var, complete_sample_ids);
			writeVariant(&output_vcf, cur_var);
		}

		if ((num_variants % 100000) == 0) {

			std::cout << "[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants" << endl;
		}

		delete cur_var;
	}

	cout << "\n[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants and " << num_alt_alleles << " alternative alleles (" << num_missing_alleles << " excluded missing alleles)." <<  endl;
	cout << "\n[" << Utils::getLocalTime() << "] Filtered " << num_filtered_alleles << " alternative alleles with an allele posterior less than " << min_allele_posterior << "." << endl;
	cout << "[" << Utils::getLocalTime() << "] Filtered " << num_filtered_samples << " samples with a genotype posterior less than " << min_genotype_posterior << "." << endl;

	cout << "\n[" << Utils::getLocalTime() << "] Wrote " << num_called_split_variants << " called bi-allelic variants." << endl;

	cout << endl;

	return 0;
}
