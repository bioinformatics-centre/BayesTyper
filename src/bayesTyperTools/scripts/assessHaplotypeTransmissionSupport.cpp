
/*
assessHaplotypeTransmissionSupport.cpp - This file is part of BayesTyper (v1.1)


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
#include <set>
#include <algorithm>
#include <iostream>
#include <fstream>

#include "VcfFile.hpp"
#include "Variant.hpp"
#include "JoiningString.hpp"
#include "Auxiliaries.hpp"

typedef vector<vector<uint> > HaplotypeConfiguration;

struct IntervalHapConf {

	uint start;
	uint end;
	HaplotypeConfiguration hap_conf;
};

typedef unordered_map<string, vector<IntervalHapConf> > HapTranIdx;

void enumerationRecursion(uint depth, uint cardinality, uint complexity, vector<uint> path, vector<vector<uint> > & all_paths) {

	uint local_depth = depth + 1;

	if (local_depth <= cardinality) {

		for (uint i = 0; i < complexity; i++) {

			vector<uint> local_path = path;
			local_path.push_back(i);
			enumerationRecursion(local_depth, cardinality, complexity, local_path, all_paths);
		}

	} else {

		all_paths.push_back(path);
	}
}

vector<vector<uint> > enumerateCombinations(uint cardinality, uint complexity) {

	vector<uint> start_path;
	vector<vector<uint> > all_paths;

	enumerationRecursion(0, cardinality, complexity, start_path, all_paths);

	return all_paths;
}

uint hapCharToIdx(char hap_char) {

	switch (hap_char) {

		case 'A' : return 0;
		case 'B' : return 1;
		case 'C' : return 2;
		case 'D' : return 3;
		default : assert(false);
	}
}

int main(int argc, char const *argv[]) {

	if (argc != 4) {

		std::cout << "USAGE: assessHaplotypeTransmissionSupport <callset.vcf> <haplotype_transmissions.txt> <output_prefix>" << std::endl;
		return 1;
	}

	cout << "\n[" << Utils::getLocalTime() << "] Running BayesTyperTools (" << BT_VERSION << ") assessHaplotypeTranmissionSupport script ...\n" << endl;

	GenotypedVcfFileReader callset_vcf_reader(argv[1], false);

	auto output_meta_data = callset_vcf_reader.metaData();
	assert(!output_meta_data.infoDescriptors().count("HTV"));

	vector<pair<string, string> > htv_descriptor_elems({make_pair("ID","HTV"), make_pair("Number","1"), make_pair("Type","String"), make_pair("Description","Variant validated by haplotype transmission (TRUE, FALSE, NA).")});
	output_meta_data.infoDescriptors().emplace("HTV", Attribute::DetailedDescriptor(htv_descriptor_elems));

	VcfFileWriter output_vcf(string(argv[3]) + ".vcf", output_meta_data, false);

	Variant * cur_callset_var;
	cout << "\n[" << Utils::getLocalTime() << "] Building index of haplotype transmissions ...\n" << endl;

	HapTranIdx haplotype_transmission_idx;
	ifstream hap_tran_file(argv[2]);
	assert(hap_tran_file.is_open());
	unordered_map<string,uint> sampleId_to_hap_idx;
	string hap_tran_line;

	uint num_parsed_intervals = 0;
	while (std::getline(hap_tran_file, hap_tran_line)) {

		auto hap_tran_line_split = Utils::splitString(hap_tran_line, '\t');

		if (hap_tran_line.front() == '#') {

			for (uint header_idx = 3; header_idx < hap_tran_line_split.size(); ++header_idx) {

				sampleId_to_hap_idx.emplace(hap_tran_line_split.at(header_idx), header_idx - 3);
			}

			assert(sampleId_to_hap_idx.size() == callset_vcf_reader.metaData().numSamples());

		} else {

			++num_parsed_intervals;

			HaplotypeConfiguration hap_conf(sampleId_to_hap_idx.size());

			for (uint sample_idx = 3; sample_idx < hap_tran_line_split.size(); ++sample_idx) {

				assert(hap_tran_line_split.at(sample_idx).size() == 1 or hap_tran_line_split.at(sample_idx).size() == 2);

				for (auto & hap_indicator : hap_tran_line_split.at(sample_idx)) {

					hap_conf.at(sample_idx-3).emplace_back(hapCharToIdx(hap_indicator));
				}
			}

			IntervalHapConf interval_hap_conf = {static_cast<uint>(stoi(hap_tran_line_split.at(1))), static_cast<uint>(stoi(hap_tran_line_split.at(2))), hap_conf};
			auto find_res = haplotype_transmission_idx.find(hap_tran_line_split.front());

			if (find_res != haplotype_transmission_idx.end()) {

				assert(find_res->second.back().end < interval_hap_conf.start);
				find_res->second.push_back(interval_hap_conf);

			} else {

				assert(haplotype_transmission_idx.emplace(hap_tran_line_split.front(), vector<IntervalHapConf>{interval_hap_conf}).second);
			}
		}
  	}

	auto tot_num_intervals = std::accumulate(haplotype_transmission_idx.begin(), haplotype_transmission_idx.end(), 0, [](int a, const HapTranIdx::value_type & b){return a + b.second.size();});
	cout << "\n[" << Utils::getLocalTime() << "] Completed haplotype transmission index with " << haplotype_transmission_idx.size() << " chromosome(s) and a total of " << tot_num_intervals << " interval(s).\n" << endl;

	cout << "\n[" << Utils::getLocalTime() << "] Assessing haplotype transmission ...\n" << endl;

	uint num_vars = 0;
	uint num_vars_outside_block = 0;
	uint num_vars_skipped_filter = 0;

	uint num_vars_pass = 0;
	uint num_vars_pass_filtered_samples = 0;
	uint num_vars_false = 0;
	uint num_vars_false_ploidy = 0;

	while (callset_vcf_reader.getNextVariant(&cur_callset_var)) {

		++num_vars;

		auto chrom_fi = haplotype_transmission_idx.find(cur_callset_var->chrom());
		if (chrom_fi != haplotype_transmission_idx.end()) {

			auto hap_block_iter = lower_bound(chrom_fi->second.begin(), chrom_fi->second.end(),
										 cur_callset_var->pos(),
										 [](const vector<IntervalHapConf>::value_type a, uint b){ return (a.end < b); });


			if (hap_block_iter != chrom_fi->second.end() and hap_block_iter->start <= cur_callset_var->pos() and (cur_callset_var->pos() + cur_callset_var->ref().seq().size() - 1) <= hap_block_iter->end) {

				assert(cur_callset_var->pos() <= hap_block_iter->end);

				vector<vector<uint> > hap_allele_assignment_combs = enumerateCombinations(4, cur_callset_var->numAlls());
				assert(!hap_allele_assignment_combs.empty());

				bool all_samples_validated;
				bool ploidy_mismatch;
				uint filtered_samples;

				for (auto & hap_allele_assignment_comb : hap_allele_assignment_combs) {

					all_samples_validated = true;
					ploidy_mismatch = false;
					filtered_samples = 0;

					for (auto & sample_id : cur_callset_var->sampleIds()) {

						if (cur_callset_var->getSample(sample_id).callStatus() == Sample::CallStatus::Complete) {

							auto sample_hap_conf = hap_block_iter->hap_conf.at(sampleId_to_hap_idx.at(sample_id));
							assert(sample_hap_conf.size() <= 2);

							vector<uint> expected_sample_gt;

							for (auto & hap_id : sample_hap_conf) {

								expected_sample_gt.emplace_back(hap_allele_assignment_comb.at(hap_id));
							}

							sort(expected_sample_gt.begin(), expected_sample_gt.end());

							auto obs_sample_gt = cur_callset_var->getSample(sample_id).genotypeEstimate();

							sort(obs_sample_gt.begin(), obs_sample_gt.end());

							if (expected_sample_gt.size() != obs_sample_gt.size()) {

								all_samples_validated = false;
								ploidy_mismatch = true;
								break;

							} else {

								if (!std::equal(expected_sample_gt.begin(), expected_sample_gt.end(), obs_sample_gt.begin())) {

									all_samples_validated = false;
									break;
								}
							}

						} else {

							++filtered_samples;
						}
					}

					if (all_samples_validated) {

						break;
					}
				}

				if (all_samples_validated) {

					assert(!ploidy_mismatch);
					cur_callset_var->info().setValue<string>("HTV", "TRUE");
					++num_vars_pass;
					if (filtered_samples > 0) {

						++num_vars_pass_filtered_samples;
					}
					// cout << "VALIDATED" << endl;

				} else {

					cur_callset_var->info().setValue<string>("HTV", "FALSE");

					if (ploidy_mismatch) {

						++num_vars_false_ploidy;
						// cout << "PLOIDY MISMATCH" << endl;
						// cout << "NEW ELIGIBLE VARIANT" << endl;
						// cout << cur_callset_var->vcf(callset_vcf_reader.metaData()) << endl;
						// cout << "Hap conf " << hap_block_iter->hap_conf << endl;

					} else {

						++num_vars_false;
						// cout << "NON-VALIDATED" << endl;
					}
				}

			} else {

				cur_callset_var->info().setValue<string>("HTV", "NA");
				++num_vars_outside_block;
			}

		} else {

			cur_callset_var->info().setValue<string>("HTV", "NA");
			++num_vars_outside_block;
		}

		output_vcf.write(cur_callset_var);
		delete cur_callset_var;
	}

	uint num_eligible_variants = num_vars - num_vars_skipped_filter - num_vars_outside_block;

	cout << "\n[" << Utils::getLocalTime() << "] Parsed " << num_vars << " variant(s)" <<  endl;
	cout << "\n[" << Utils::getLocalTime() << "] " << num_vars_outside_block << " variant(s) were outside a haplotype block" <<  endl;
	cout << "\n[" << Utils::getLocalTime() << "] leaving " <<  num_eligible_variants << " variant(s) eligible for haplotype transmission evaluation" <<  endl;

	cout << "\n[" << Utils::getLocalTime() << "] " << num_vars_pass << " eligible variant(s) were supported by haplotype transmission" <<  endl;
	cout << "\n[" << Utils::getLocalTime() << "] of these " << num_vars_pass_filtered_samples << " eligible variant(s) had filtered samples " <<  endl;
	cout << "\n[" << Utils::getLocalTime() << "] " << num_vars_false + num_vars_false_ploidy << " eligble variant(s) were NOT supported by haplotype transmission" <<  endl;
	cout << "\n[" << Utils::getLocalTime() << "] of these " << num_vars_false_ploidy << " eligible variant(s) had were NOT supported due to ploidy mismatch" <<  endl;

	return 0;
}
