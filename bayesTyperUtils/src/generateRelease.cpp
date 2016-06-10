
/*
generateRelease.cpp - This file is part of BayesTyper (v0.9)


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


#include <unordered_map>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <math.h>
#include <unordered_set>
#include <algorithm>

#include "vcf++/JoiningString.hpp"
#include "vcf++/Auxiliaries.hpp"
#include "vcf++/Stats.hpp"
#include "vcf++/JoiningString.hpp"

#include "generateRelease.hpp"


namespace GenerateRelease {

	void generateRelease(const string & vcf_filename, const string & output_prefix, const float min_called_probability) {

		cout << "[" << Utils::getLocalTime() << "] Running BayesTyperUtils (" << BTU_VERSION << ") generateRelease ...\n" << endl;

		GenotypedVcfFileReader vcf_reader(vcf_filename, true);

		assert(vcf_reader.metaData().miscMeta().erase("CommandLine"));

		assert(vcf_reader.metaData().infoDescriptors().erase("AE"));
		assert(vcf_reader.metaData().infoDescriptors().erase("ANC"));
		assert(vcf_reader.metaData().infoDescriptors().erase("HCR"));
		assert(vcf_reader.metaData().infoDescriptors().erase("HRS"));
		assert(vcf_reader.metaData().infoDescriptors().erase("VCGI"));
		assert(vcf_reader.metaData().infoDescriptors().erase("VCGS"));
		assert(vcf_reader.metaData().infoDescriptors().erase("VCI"));
		assert(vcf_reader.metaData().infoDescriptors().erase("VCS"));
		assert(vcf_reader.metaData().infoDescriptors().erase("VT"));

		auto complete_meta_data = vcf_reader.metaData();

		assert(complete_meta_data.infoDescriptors().emplace("AT", Attribute::DetailedDescriptor("AT", Attribute::Number::A, Attribute::Type::String, "Allele type")).second);

		string lvq_label("LVQ");
		assert(complete_meta_data.filterDescriptors().emplace(lvq_label, Attribute::Descriptor(lvq_label, "All alternative alleles excluding missing have low allele call probability (ACP < " + to_string(min_called_probability) + ")")).second);

		string nc_label("NC");
		assert(complete_meta_data.filterDescriptors().emplace(nc_label, Attribute::Descriptor(nc_label, "All alternative alleles excluding missing are not called (ACP == 0)")).second);

		VcfMetaData called_file_meta_data = complete_meta_data;
		called_file_meta_data.filterDescriptors().clear();
		called_file_meta_data.clearSamples();

		assert(called_file_meta_data.infoDescriptors().erase("ACP"));				
		assert(called_file_meta_data.infoDescriptors().erase("ACO"));				
		
		called_file_meta_data.infoDescriptors().erase("AsmVar_ASQR");
		called_file_meta_data.infoDescriptors().erase("DNE");
		called_file_meta_data.infoDescriptors().erase("RMA");

		VcfMetaData de_novo_file_meta_data = complete_meta_data;
		de_novo_file_meta_data.filterDescriptors().clear();

		VcfFileWriter vcf_writer_complete(output_prefix + "_complete.vcf", complete_meta_data, true);
		VcfFileWriter vcf_writer_called(output_prefix + "_called_alleles.vcf", called_file_meta_data, true);
		VcfFileWriter vcf_writer_de_novo(output_prefix + "_de_novo.vcf", de_novo_file_meta_data, true);

		Variant * cur_var;

		uint num_variants = 0;
		uint num_alleles = 0;

		uint num_low_prob_variants = 0;
		uint num_filt_variants = 0;

		uint num_missing_alleles = 0;
		uint num_not_called_variants = 0;
		uint num_low_prob_alleles = 0;

		uint num_de_novo_events = 0;
		uint num_de_novo_multiple_events = 0;

		assert(min_called_probability > 0);
		assert(min_called_probability <= 1);

		while (vcf_reader.getNextVariant(&cur_var)) {

			num_variants++;
			num_alleles += cur_var->numAlls();

			assert(cur_var->filters().size() == 1);

			// Label alleles
			vector<uint> alt_rm_idxs;
			alt_rm_idxs.reserve(cur_var->numAlts());

			float alt_acp_max = 0;

			for (uint allele_idx = 0; allele_idx < cur_var->numAlls(); allele_idx++) {

				assert(cur_var->info().hasValue("AN"));

				auto acp_value = cur_var->allele(allele_idx).info().getValue<float>("ACP");
				assert(acp_value.second);

				if (allele_idx > 0) {

					assert(cur_var->allele(allele_idx).info().hasValue("AC"));
					assert(cur_var->allele(allele_idx).info().hasValue("AF"));
					assert(cur_var->allele(allele_idx).info().hasValue("ACO"));
					assert(cur_var->allele(allele_idx).info().hasValue("AAI"));

					auto at = Auxiliaries::alleleAttributes(cur_var->allele(allele_idx), cur_var->ref());
					auto rma_value = cur_var->allele(allele_idx).info().getValue<string>("RMA");

					if (rma_value.second) {

						assert((rma_value.first == ".") or (at.type == Auxiliaries::Type::Insertion) or (at.type == Auxiliaries::Type::Deletion));
					}

					assert(cur_var->allele(allele_idx).info().setValue("AT", at.typeStr()));

					if (cur_var->allele(allele_idx).isMissing()) {

						alt_rm_idxs.push_back(allele_idx - 1);
						num_missing_alleles++;

						continue;
					}

					alt_acp_max = max(alt_acp_max, acp_value.first);
					
					if (acp_value.first < min_called_probability) {

						alt_rm_idxs.push_back(allele_idx - 1);
						num_low_prob_alleles++;
					}
				}
			}

			if (cur_var->filters().front() == "PASS") {

				auto var_qual = cur_var->qual();
				assert(var_qual.second);

				auto var_call_prob = 1 - pow(10, -1 * var_qual.first/10);
				assert((var_call_prob == alt_acp_max) or (abs(var_call_prob - alt_acp_max) < pow(10, -3)));

				if (floatCompare(alt_acp_max, 0)) {

					assert(alt_rm_idxs.size() == cur_var->numAlts());

					cur_var->setFilters({nc_label});
					num_not_called_variants++;
			
				} else if (alt_acp_max < min_called_probability) {

					assert(alt_rm_idxs.size() == cur_var->numAlts());

					cur_var->setFilters({lvq_label});
					num_low_prob_variants++;

				} else {

					assert(alt_rm_idxs.size() < cur_var->numAlts());
				}

			} else {

				num_filt_variants++;
			}

			vcf_writer_complete.write(cur_var);

			auto dne_att = cur_var->info().getValue<string>("DNE");

			if (dne_att.second) {

				assert(cur_var->filters().front() == "PASS");

				if (dne_att.first.find(",") == string::npos) {

					num_de_novo_events++;

				} else {

					num_de_novo_multiple_events++;
				}

				vcf_writer_de_novo.write(cur_var);
			}

			if (cur_var->filters().front() == "PASS") {

				assert(alt_rm_idxs.size() < cur_var->numAlts());

				int alt_rm_ac = 0;

				for (auto & alt_rm_idx: alt_rm_idxs) {

					auto ac_value = cur_var->alt(alt_rm_idx).info().getValue<int>("AC");
					assert(ac_value.second);

					alt_rm_ac += ac_value.first;
				}

				cur_var->clearSamples();
				cur_var->removeAlts(alt_rm_idxs);
	            
	            auto an_value = cur_var->info().getValue<int>("AN");
	            assert(an_value.second);
	            assert(an_value.first >= alt_rm_ac);

	            if (alt_rm_ac > 0) {

		            assert(!cur_var->info().setValue<int>("AN", an_value.first - alt_rm_ac));
					assert(cur_var->numAlts() > 0);

					for (uint alt_idx = 0; alt_idx < cur_var->numAlts(); alt_idx++) {

		                if (an_value.first == alt_rm_ac) {

			                assert(!cur_var->alt(alt_idx).info().setValue<float>("AF", 0));

		                } else {

		                	auto ac_value = cur_var->alt(alt_idx).info().getValue<int>("AC");
		                	assert(ac_value.second);
			                
			                assert(!cur_var->alt(alt_idx).info().setValue<float>("AF", ac_value.first / static_cast<float>(an_value.first - alt_rm_ac)));
		                }
		            }
	        	}

				vcf_writer_called.write(cur_var);
			}

			delete cur_var;

			if ((num_variants % 100000) == 0) {

				std::cout << "[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants" << endl;
			}
		}

		cout << "\n[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants and " << num_alleles << " alleles (including reference and missing)\n" << endl;
		cout << "\t\t - " << num_filt_variants << " variant(s) were already labelled as filtered (UV or ISD)" << endl;
		cout << "\t\t - " << num_not_called_variants << " variant(s) were labelled as not called (NC)" << endl;
		cout << "\t\t - " << num_low_prob_variants << " variant(s) were labelled as low variant quality (LVQ)\n" << endl;
		cout << "\t\t - " << num_missing_alleles << " missing alleles(s)" << endl;
		cout << "\t\t - " << num_low_prob_alleles << " alleles(s) had probability below threshold\n" << endl;
		cout << "\t\t - " << num_de_novo_events << " de novo events" << endl;
		cout << "\t\t - " << num_de_novo_multiple_events << " multi de novo events" << endl;
		cout << endl;
	}
}
