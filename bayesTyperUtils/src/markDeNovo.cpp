
/*
markDeNovo.cpp - This file is part of BayesTyper (v0.9)


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


#include <memory>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <algorithm>

#include "vcf++/JoiningString.hpp"
#include "vcf++/Trio.hpp"
#include "vcf++/Auxiliaries.hpp"

#include "markDeNovo.hpp"

namespace MarkDeNovo {

	static const vector<string> de_novo_exclusion_attributes({"HCR","AE","ANC"});

	void updatePositiveAllelePosteriorCount(vector<uint> * positive_allele_posterior_count, Sample & sample) {

		assert(positive_allele_posterior_count->size() == sample.alleleInfo().size());

		if (sample.isInformative()) {

			assert(sample.callStatus() == Sample::CallStatus::Complete);
			assert(sample.ploidy() == Sample::Ploidy::Diploid);

			for (uint i = 0; i < positive_allele_posterior_count->size(); i++) {

				auto map_value = sample.alleleInfo().at(i).getValue<float>("MAP");
				assert(map_value.second);
				
				if (!(Utils::floatCompare(map_value.first, 0))) {

					positive_allele_posterior_count->at(i)++;
				}	
			}

		} else {

			assert(sample.callStatus() == Sample::CallStatus::Missing);

			for (uint i = 0; i < positive_allele_posterior_count->size(); i++) {

				positive_allele_posterior_count->at(i)++;
			}
		}
	}

	void markDeNovo(const string & vcf_filename, const string & output_prefix, const string & trio_info_str, const bool use_genome_dk_trio_syntax, const float de_novo_zero_inflation_threshold, const float de_novo_allelic_balance_deviation) {

        cout << "[" << Utils::getLocalTime() << "] Running BayesTyperUtils (" << BTU_VERSION << ") markDeNovo ...\n" << endl;

		GenotypedVcfFileReader vcf_reader(vcf_filename, true);		
	    	
    	vector<pair<string, string> > dne_descriptor_elems({make_pair("ID","DNE"), make_pair("Number","."), make_pair("Type","String"), make_pair("Description","De novo event(s) (<child_id>:<ancestral_allele_idx>:<de_novo_allele_idx>:<event_type>:<child_allelic_balance>,...).")});
		assert(vcf_reader.metaData().infoDescriptors().emplace("DNE", Attribute::DetailedDescriptor(dne_descriptor_elems)).second);

		VcfFileWriter vcf_writer_dne(output_prefix + ".vcf", vcf_reader.metaData(), true);
	
		vector<Trio::TrioInfo> all_trio_info;

		if (use_genome_dk_trio_syntax) {

			assert(trio_info_str.empty());
			all_trio_info = Trio::parseGenomeDKPedigree(vcf_reader.metaData());

		} else {

			all_trio_info = Trio::parsePedigree(vcf_reader.metaData(), trio_info_str);
		}

		cout << "[" << Utils::getLocalTime() << "] Marking de novo events between:\n" << endl;

		for (auto &cur_trio_info: all_trio_info) {

			cout << "\t- " << cur_trio_info.id << ": Father " << cur_trio_info.father << ", mother " << cur_trio_info.mother << ", and child " << cur_trio_info.child << endl;
		}

		cout << endl;

		uint num_variants = 0;
		uint num_filtered = 0;
		uint num_non_autosomal = 0;
		uint num_filtered_trios = 0;

		uint num_trio_bi = 0;
		uint num_trio_multi = 0;

		uint num_de_novo_events_bi = 0;
		uint num_de_novo_events_multi = 0;

		unordered_map<string, uint> trio_coverage_stats;

		Variant * cur_var;

		while (vcf_reader.getNextVariant(&cur_var)) {

			num_variants++;

			assert(cur_var->filters().size() == 1);

			if (cur_var->filters().front() != "PASS") {

				num_filtered++;
			
			} else if (vcf_reader.metaData().getContig(cur_var->chrom()).type() != Contig::Type::Autosomal) {

				num_non_autosomal++;
			
			} else {

				string variant_type = Auxiliaries::variantType(*cur_var);
					
				uint num_non_zero_prob_alleles = Auxiliaries::getNonZeroProbAlleleIdxsSorted(*cur_var).size();
				assert(num_non_zero_prob_alleles <= cur_var->numAlls());

				vector<Trio> de_novo_trios;
				vector<pair<string, Trio> > coverage_stats_trios;
				coverage_stats_trios.reserve(all_trio_info.size());

				vector<uint> positive_allele_posterior_count(cur_var->numAlls(), 0);

				for (auto &sample_id: cur_var->sampleIds()) {

					updatePositiveAllelePosteriorCount(&positive_allele_posterior_count, cur_var->getSample(sample_id));
				}

				for (auto &cur_trio_info: all_trio_info) {

					Trio cur_trio = Trio(*cur_var, cur_trio_info);
					assert(cur_trio.isDiploid());

					if (cur_trio.isFiltered()) {

						num_filtered_trios++;

					} else {

						if (cur_var->numAlls() ==  2) {

							num_trio_bi++;
						
						} else {

							num_trio_multi++;
						}

						bool is_proper_event = cur_trio.isExclusivelyChildHeterozygote();

						if (!(Utils::floatCompare(cur_trio.minGPP(), 1))) {

							is_proper_event = false;
						} 

						for (auto &att: de_novo_exclusion_attributes) {

							if (cur_var->info().hasValue(att)) {

								is_proper_event = false;
								break;
							}
						}

						if (Auxiliaries::hasAmbiguous(*cur_var)) {

							is_proper_event = false;
						}

						if (cur_trio.isConcordant()) {

							if (is_proper_event and !cur_trio.hasCalledMissing()) {

								coverage_stats_trios.emplace_back(cur_trio_info.id, cur_trio);
							}					

						} else if (cur_trio.isDeNovo()) {

							auto de_novo_event = cur_trio.deNovoEvent();					

							auto aai_att = cur_var->allele(de_novo_event.de_novo_allele_idx).info().getValue<string>("AAI"); 
							assert(aai_att.second);

							if (aai_att.first != ".") {

								auto aai_att_split = Utils::splitString(aai_att.first, ':');
								
								assert(!(aai_att_split.empty()));
								assert(find(aai_att_split.begin(), aai_att_split.end(), ".") == aai_att_split.end());
							
							} else if (is_proper_event) {

								assert(positive_allele_posterior_count.at(de_novo_event.de_novo_allele_idx) >= 1);

								if (positive_allele_posterior_count.at(de_novo_event.de_novo_allele_idx) == 1) {

									if (Auxiliaries::hasMissing(*cur_var)) {

										for (uint all_idx = 0; all_idx < (cur_var->numAlls() - 1); all_idx++) {

											assert(!(cur_var->allele(all_idx).isMissing()));
										}

										assert(cur_var->allele(cur_var->numAlls() - 1).isMissing());

										if (positive_allele_posterior_count.back() == 0) {

											coverage_stats_trios.emplace_back(cur_trio_info.id, cur_trio);
											de_novo_trios.emplace_back(cur_trio);	
										} 
									
									} else {

										coverage_stats_trios.emplace_back(cur_trio_info.id, cur_trio);
										de_novo_trios.emplace_back(cur_trio);
									}
								}
							}
						}
					}
				}

				if (!(de_novo_trios.empty())) {

					JoiningString dne_elements(',');

					for (auto &de_novo_trio: de_novo_trios) {

						assert(de_novo_trio.isDeNovo());
						auto de_novo_event = de_novo_trio.deNovoEvent();					

						assert(de_novo_event.de_novo_allele_idx > 0);
						assert(positive_allele_posterior_count.at(de_novo_event.de_novo_allele_idx) == 1);

						float max_zero_inflation = 1;
						float child_allelic_balance = de_novo_trio.childAllelicBalance();

						if (!(Utils::floatCompare(max_zero_inflation, de_novo_zero_inflation_threshold)) and (max_zero_inflation > de_novo_zero_inflation_threshold)) {

							continue;							
						}

						float cur_deviation = abs(child_allelic_balance - 0.5);

						if (!(Utils::floatCompare(cur_deviation, de_novo_allelic_balance_deviation)) and (cur_deviation > de_novo_allelic_balance_deviation)) {

							continue;
						}

						assert(num_non_zero_prob_alleles > 1);
				
						if (num_non_zero_prob_alleles ==  2) {

							num_de_novo_events_bi++;
						
						} else {

							num_de_novo_events_multi++;
						}

						dne_elements.join(de_novo_event.child_id + ":" + to_string(de_novo_event.ancestral_allele_idx) + ":" + to_string(de_novo_event.de_novo_allele_idx) + ":" + Auxiliaries::alleleAttributes(cur_var->allele(de_novo_event.de_novo_allele_idx), cur_var->allele(de_novo_event.ancestral_allele_idx)).typeStr() + ":" + to_string(child_allelic_balance));
					}

					if (!(dne_elements.empty())) {

						assert(cur_var->info().setValue<string>("DNE", dne_elements.str()));
					}
				}

				for (auto &stats_trio: coverage_stats_trios) {
					
					assert(num_non_zero_prob_alleles > 1);

					JoiningString trio_coverage_stats_line_elements('\t');
					trio_coverage_stats_line_elements.join(stats_trio.first);		
					trio_coverage_stats_line_elements.join(Utils::boolToString(stats_trio.second.isDeNovo()));		
					trio_coverage_stats_line_elements.join(Utils::floatToString(stats_trio.second.minNumKmers(), 0));
					trio_coverage_stats_line_elements.join(Utils::floatToString(stats_trio.second.minNumObservedKmers(), 0));
					trio_coverage_stats_line_elements.join(Utils::floatToString(stats_trio.second.minNumUniqueKmers(), 0));
					trio_coverage_stats_line_elements.join(Utils::floatToString(stats_trio.second.minObservedKmerCoverage(), 1));
					trio_coverage_stats_line_elements.join(Utils::floatToString(stats_trio.second.childAllelicBalance(), 2));
					trio_coverage_stats_line_elements.join(variant_type);
					trio_coverage_stats_line_elements.join(to_string(num_non_zero_prob_alleles));
	
					auto trio_coverage_stats_emplace = trio_coverage_stats.emplace(trio_coverage_stats_line_elements.str(), 0);
					trio_coverage_stats_emplace.first->second++;
				}
			}

			vcf_writer_dne.write(cur_var);
			
			delete cur_var;

			if ((num_variants % 100000) == 0) {

				std::cout << "[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants" << endl;
			}
		}

	    ofstream trio_coverage_stats_writer(output_prefix + "_coverage_stats.txt");
	    assert(trio_coverage_stats_writer.is_open());

		trio_coverage_stats_writer << "Count\tTrioId\tIsDeNovo\tMinTrioNumKmers\tMinTrioNumObservedKmers\tMinTrioNumUniqueKmers\tMinTrioObservedKmerCoverage\tChildAllelicBalance\tVariantType\tNumNonZeroProbAlleles\n";

		for (auto &stats: trio_coverage_stats) {

			trio_coverage_stats_writer << stats.second << "\t" << stats.first << "\n";
		}

		trio_coverage_stats_writer.close();

		cout << "\n[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants:\n" << endl;
		cout << "\t\t- Skipped " << num_filtered << " filtered variants" << endl;
		cout << "\t\t- Skipped " << num_non_autosomal << " variants from non autosomal chromosomes" << endl;

		uint remaining_variants = num_variants - num_filtered - num_non_autosomal;

		cout << "\n[" << Utils::getLocalTime() << "] Parsed " << remaining_variants * all_trio_info.size() << " trio genotypes:\n" << endl;
		cout << "\t\t- Skipped " << num_filtered_trios << " filtered trio genotypes (at least one genotype filtered)" << endl;

		cout << "\n[" << Utils::getLocalTime() << "] Of the remaining " << num_trio_bi + num_trio_multi << " trio genotypes " << num_de_novo_events_bi + num_de_novo_events_multi << " had a de novo event:\n"<< endl;
		cout << "\t\t- Bi-allelic: " << num_de_novo_events_bi << " (" << static_cast<float>(num_de_novo_events_bi)/num_trio_bi << ")" << endl;
		cout << "\t\t- Multi-allelic: " << num_de_novo_events_multi << " (" << static_cast<float>(num_de_novo_events_multi)/num_trio_multi << ")" << endl;

		cout << "\n" << endl;
	}
}
