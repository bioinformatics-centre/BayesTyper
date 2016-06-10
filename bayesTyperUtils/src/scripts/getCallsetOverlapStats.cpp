
/*
getCallsetOverlapStats.cpp - This file is part of BayesTyper (v0.9)


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
#include <algorithm>
#include <unordered_set>

#include "vcf++/VcfFile.hpp"
#include "vcf++/Variant.hpp"
#include "vcf++/JoiningString.hpp"
#include "vcf++/Trio.hpp"
#include "vcf++/Auxiliaries.hpp"


string getFilterString(Variant * cur_var) {

	if (!cur_var) {

		return ".";
	}

	if (cur_var->filters().empty()) {

		return ".";

	} else {

		JoiningString filter_elements(';');
		filter_elements.join(cur_var->filters());

		return filter_elements.str();
	}
}

string getCallProbabilityString(Variant * cur_var, const uint alt_idx) {

	if (!cur_var) {

		return "0";
	}

    auto acp_value = cur_var->alt(alt_idx).info().getValue<float>("ACP");

	if (acp_value.second) {
	
		return Utils::floatToString(acp_value.first, 2);

	} else {

		return "1";
	}
}

void addOverlapStats(unordered_map<string, ulong> * allele_overlap_stats, const string & chromosome_type, Variant * cur_var_1, const uint alt_idx_1, Variant * cur_var_2, const uint alt_idx_2) {

	assert(cur_var_1 or cur_var_2);

	JoiningString allele_overlap_stats_line_elements('\t');
	allele_overlap_stats_line_elements.join(chromosome_type);

	allele_overlap_stats_line_elements.join(getFilterString(cur_var_1));
	allele_overlap_stats_line_elements.join(getCallProbabilityString(cur_var_1, alt_idx_1));

	allele_overlap_stats_line_elements.join(getFilterString(cur_var_2));
	allele_overlap_stats_line_elements.join(getCallProbabilityString(cur_var_2, alt_idx_2));

	if (cur_var_1 and cur_var_2) {

		auto allele_attributes_1 = Auxiliaries::alleleAttributes(cur_var_1->alt(alt_idx_1), cur_var_1->ref());
		auto allele_attributes_2 = Auxiliaries::alleleAttributes(cur_var_2->alt(alt_idx_2), cur_var_2->ref());

		assert(allele_attributes_1.typeStr() == allele_attributes_2.typeStr());
		assert(allele_attributes_1.length == allele_attributes_2.length);
		assert(allele_attributes_1.num_ambiguous == allele_attributes_2.num_ambiguous);
		assert(allele_attributes_1.sv_length == allele_attributes_2.sv_length);

		allele_overlap_stats_line_elements.join(allele_attributes_1.typeStr());
		allele_overlap_stats_line_elements.join(to_string(allele_attributes_1.length));
		allele_overlap_stats_line_elements.join(to_string(allele_attributes_1.num_ambiguous));
		allele_overlap_stats_line_elements.join(to_string(allele_attributes_1.sv_length));	

	} else if (cur_var_1) {
	
		auto allele_attributes = Auxiliaries::alleleAttributes(cur_var_1->alt(alt_idx_1), cur_var_1->ref());

		allele_overlap_stats_line_elements.join(allele_attributes.typeStr());
		allele_overlap_stats_line_elements.join(to_string(allele_attributes.length));
		allele_overlap_stats_line_elements.join(to_string(allele_attributes.num_ambiguous));
		allele_overlap_stats_line_elements.join(to_string(allele_attributes.sv_length));	

	} else {

		assert(cur_var_2);

		auto allele_attributes = Auxiliaries::alleleAttributes(cur_var_2->alt(alt_idx_2), cur_var_2->ref());

		allele_overlap_stats_line_elements.join(allele_attributes.typeStr());
		allele_overlap_stats_line_elements.join(to_string(allele_attributes.length));
		allele_overlap_stats_line_elements.join(to_string(allele_attributes.num_ambiguous));
		allele_overlap_stats_line_elements.join(to_string(allele_attributes.sv_length));
	}	

	auto allele_overlap_stats_emplace = allele_overlap_stats->emplace(allele_overlap_stats_line_elements.str(), 0);
	allele_overlap_stats_emplace.first->second++;
}

void addTrioStats(unordered_map<string, ulong> * concordance_trio_stats, const string & chromosome_type, Variant * cur_var_1, Variant * cur_var_2, const vector<Trio::TrioInfo> & all_trio_info) {

	assert(cur_var_1);
	assert(cur_var_2);

	assert(cur_var_1->chrom() == cur_var_2->chrom());
	assert(cur_var_1->pos() == cur_var_2->pos());

	const string variant_type = Auxiliaries::variantType(*cur_var_1);
	assert(variant_type == Auxiliaries::variantType(*cur_var_2));

	const bool has_missing = Auxiliaries::hasMissing(*cur_var_1);
	assert(has_missing == Auxiliaries::hasMissing(*cur_var_2));

	const bool has_repeat = (Auxiliaries::hasRepeat(*cur_var_1) or Auxiliaries::hasRepeat(*cur_var_2));

	uint max_allele_length_1 = 0;
    uint max_num_ambiguous_1 = 0;
	uint max_abs_sv_length_1 = 0;

	uint max_allele_length_2 = 0;
    uint max_num_ambiguous_2 = 0;
	uint max_abs_sv_length_2 = 0;

	assert(cur_var_1->numAlts() == cur_var_2->numAlts());
	assert(cur_var_1->numAlls() == cur_var_2->numAlls());

	for (uint alt_idx = 0; alt_idx < cur_var_1->numAlts(); alt_idx++) {

		auto allele_attributes_1 = Auxiliaries::alleleAttributes(cur_var_1->alt(alt_idx), cur_var_1->ref());
		auto allele_attributes_2 = Auxiliaries::alleleAttributes(cur_var_2->alt(alt_idx), cur_var_2->ref());
	
		max_allele_length_1 = max(max_allele_length_1, allele_attributes_1.length);
		max_num_ambiguous_1 = max(max_num_ambiguous_1, allele_attributes_1.num_ambiguous);
		max_abs_sv_length_1 = max(max_abs_sv_length_1, static_cast<uint>(abs(allele_attributes_1.sv_length)));

		max_allele_length_2 = max(max_allele_length_2, allele_attributes_2.length);
		max_num_ambiguous_2 = max(max_num_ambiguous_2, allele_attributes_2.num_ambiguous);
		max_abs_sv_length_2 = max(max_abs_sv_length_2, static_cast<uint>(abs(allele_attributes_2.sv_length)));
	}

	assert(max_allele_length_1 == max_allele_length_2);
	assert(max_num_ambiguous_1 == max_num_ambiguous_2);
	assert(max_abs_sv_length_1 == max_abs_sv_length_2);

	for (auto &cur_trio_info: all_trio_info) {

		Trio cur_trio_1 = Trio(*cur_var_1, cur_trio_info);
		Trio cur_trio_2 = Trio(*cur_var_2, cur_trio_info);

		if (cur_trio_1.isFiltered() or cur_trio_2.isFiltered()) {

			continue;
		}

		if (cur_trio_1.isInformative() and cur_trio_2.isInformative()) {

			JoiningString trio_stats_line_elements('\t');
			trio_stats_line_elements.join(chromosome_type);
			trio_stats_line_elements.join(cur_trio_info.id);
			trio_stats_line_elements.join(Utils::boolToString(cur_trio_1.isConcordant()));			
			trio_stats_line_elements.join(Utils::boolToString(cur_trio_2.isConcordant()));
			trio_stats_line_elements.join(getFilterString(cur_var_1));			
			trio_stats_line_elements.join(getFilterString(cur_var_2));
			trio_stats_line_elements.join(Utils::boolToString(cur_trio_1.isReferenceCall()));
			trio_stats_line_elements.join(Utils::boolToString(cur_trio_2.isReferenceCall()));
			trio_stats_line_elements.join(Utils::boolToString(cur_trio_1.isParentsBiAllelicHeterozygote()));
			trio_stats_line_elements.join(Utils::boolToString(cur_trio_2.isParentsBiAllelicHeterozygote()));
			trio_stats_line_elements.join(Utils::floatToString(cur_trio_1.minGPP(), 2));
			trio_stats_line_elements.join(Utils::floatToString(cur_trio_2.minGPP(), 2));
			trio_stats_line_elements.join(variant_type);
			trio_stats_line_elements.join(Utils::boolToString(has_missing));
			trio_stats_line_elements.join(Utils::boolToString(has_repeat));
			trio_stats_line_elements.join(to_string(cur_var_1->numAlls()));
			trio_stats_line_elements.join(to_string(max_allele_length_1));
			trio_stats_line_elements.join(to_string(max_num_ambiguous_1));
			trio_stats_line_elements.join(to_string(max_abs_sv_length_1));

			auto concordance_trio_stats_emplace = concordance_trio_stats->emplace(trio_stats_line_elements.str(), 0);
			concordance_trio_stats_emplace.first->second++;
		}
	}
}


int main(int argc, char const *argv[]) {

    if ((argc != 4) and (argc != 5)) {

        std::cout << "USAGE: getCallsetOverlapStats <callset1> <callset2> <output-prefix> (<trio-info (see \"calcTrioConcordance\" options for more details)>)" << std::endl;
        return 1;
    }

    cout << "\n[" << Utils::getLocalTime() << "] Running BayesTyperUtils (" << BTU_VERSION << ") getCallsetOverlapStats script ...\n" << endl;

	GenotypedVcfFileReader vcf_reader_1(argv[1], true);
	GenotypedVcfFileReader vcf_reader_2(argv[2], true);

	Auxiliaries::removeNonRelevantInfoDescriptors(&(vcf_reader_1.metaData()), {"ACP"});
	Auxiliaries::removeNonRelevantFormatDescriptors(&(vcf_reader_1.metaData()), {"GT", "GPP"});

	Auxiliaries::removeNonRelevantInfoDescriptors(&(vcf_reader_2.metaData()), {"ACP"});
	Auxiliaries::removeNonRelevantFormatDescriptors(&(vcf_reader_2.metaData()), {"GT", "GPP"});
	
	auto contigs = vcf_reader_1.metaData().contigs();
	assert(contigs == vcf_reader_2.metaData().contigs());

	vector<Trio::TrioInfo> all_trio_info;

	assert(vcf_reader_1.metaData().sampleIds() == vcf_reader_2.metaData().sampleIds());

	if (argc == 4) {

		all_trio_info = Trio::parseGenomeDKPedigree(vcf_reader_1.metaData());

	} else {

		all_trio_info = Trio::parsePedigree(vcf_reader_1.metaData(), string(argv[4]));
	}

	cout << "[" << Utils::getLocalTime() << "] Assessing concordance between:\n" << endl;

	for (auto &cur_trio_info: all_trio_info) {

		cout << "\t- " << cur_trio_info.id << ": Father " << cur_trio_info.father << ", mother " << cur_trio_info.mother << ", and child " << cur_trio_info.child << endl;
	}

	cout << endl;

	Variant * cur_var_1;
	Variant * cur_var_2;

	bool cur_var_1_parsed = vcf_reader_1.getNextVariant(&cur_var_1);
	bool cur_var_2_parsed = vcf_reader_2.getNextVariant(&cur_var_2);

	unordered_map<string, ulong> allele_overlap_stats;	
	unordered_map<string, ulong> concordance_trio_stats;

	ulong num_variants_1 = 0;
	ulong num_variants_2 = 0;
	ulong num_variants_overlap = 0;

	for (auto &contig: contigs) {

		string chromosome_type = vcf_reader_1.metaData().getContig(contig.id()).typeStr();
		assert(chromosome_type == vcf_reader_2.metaData().getContig(contig.id()).typeStr());

		while (cur_var_1_parsed) {

			if (cur_var_1->chrom() != contig.id()) {

				break;
			}

			while (cur_var_2_parsed) {

				if (cur_var_2->chrom() != contig.id()) {

					break;
				}	

				assert(cur_var_1->chrom() == cur_var_2->chrom());

				if (cur_var_1->pos() > cur_var_2->pos()) {

					for (uint alt_idx = 0; alt_idx < cur_var_2->numAlts(); alt_idx++) {

						addOverlapStats(&allele_overlap_stats, chromosome_type, nullptr, 0, cur_var_2, alt_idx);
					}

					num_variants_2++;
					delete cur_var_2;

					cur_var_2_parsed = vcf_reader_2.getNextVariant(&cur_var_2);

				} else if (cur_var_1->pos() == cur_var_2->pos()) {

					vector<pair<string, string> > allele_pairs_1;
					vector<pair<string, string> > allele_pairs_2;

					allele_pairs_1.reserve(cur_var_1->numAlts());
					allele_pairs_2.reserve(cur_var_2->numAlts());

					for (uint alt_idx = 0; alt_idx < cur_var_1->numAlts(); alt_idx++) {

		                Allele cur_ref_allele_1 = cur_var_1->ref();
		                Allele cur_alt_allele_1 = cur_var_1->alt(alt_idx);

		                Auxiliaries::rightTrimAllelePair(&cur_ref_allele_1, &cur_alt_allele_1);
		                allele_pairs_1.emplace_back(cur_ref_allele_1.seq(), cur_alt_allele_1.seq());
		            }

					for (uint alt_idx = 0; alt_idx < cur_var_2->numAlts(); alt_idx++) {

		                Allele cur_ref_allele_2 = cur_var_2->ref();
		                Allele cur_alt_allele_2 = cur_var_2->alt(alt_idx);

		                Auxiliaries::rightTrimAllelePair(&cur_ref_allele_2, &cur_alt_allele_2);
		                allele_pairs_2.emplace_back(cur_ref_allele_2.seq(), cur_alt_allele_2.seq());
		            }

		            assert(allele_pairs_1.size() == cur_var_1->numAlts());
		            assert(allele_pairs_2.size() == cur_var_2->numAlts());

					unordered_set<uint> has_overlap_alt_idx_2;

					for (uint alt_idx_1 = 0; alt_idx_1 < cur_var_1->numAlts(); alt_idx_1++) {

						bool has_overlap = false;

						for (uint alt_idx_2 = 0; alt_idx_2 < cur_var_2->numAlts(); alt_idx_2++) {

							if (allele_pairs_1.at(alt_idx_1) == allele_pairs_2.at(alt_idx_2)) {

								addOverlapStats(&allele_overlap_stats, chromosome_type, cur_var_1, alt_idx_1, cur_var_2, alt_idx_2);

								assert(has_overlap_alt_idx_2.insert(alt_idx_2).second);
								has_overlap = true;

								break;
							}
						}

						if (!has_overlap) {

							addOverlapStats(&allele_overlap_stats, chromosome_type, cur_var_1, alt_idx_1, nullptr, 0);
						}
					}

					for (uint alt_idx_2 = 0; alt_idx_2 < cur_var_2->numAlts(); alt_idx_2++) {

						if (has_overlap_alt_idx_2.count(alt_idx_2) > 0) {

							continue;
						}

						addOverlapStats(&allele_overlap_stats, chromosome_type, nullptr, 0, cur_var_2, alt_idx_2);
					}

		            sort(allele_pairs_1.begin(), allele_pairs_1.end());
		            sort(allele_pairs_2.begin(), allele_pairs_2.end());

		            assert(allele_pairs_1.size() == cur_var_1->numAlts());
		            assert(allele_pairs_2.size() == cur_var_2->numAlts());

					if (allele_pairs_1 == allele_pairs_2) {

						num_variants_overlap++;						
						addTrioStats(&concordance_trio_stats, chromosome_type, cur_var_1, cur_var_2, all_trio_info);
					}

					num_variants_2++;
					delete cur_var_2;

					cur_var_2_parsed = vcf_reader_2.getNextVariant(&cur_var_2);
					break;

				} else {

					for (uint alt_idx = 0; alt_idx < cur_var_1->numAlts(); alt_idx++) {

						addOverlapStats(&allele_overlap_stats, chromosome_type, cur_var_1, alt_idx, nullptr, 0);
					}

					break;
				}
			}

			num_variants_1++;
			delete cur_var_1;

			cur_var_1_parsed = vcf_reader_1.getNextVariant(&cur_var_1);
		}

		while (cur_var_2_parsed) {

			if (cur_var_2->chrom() != contig.id()) {

				break;
			}

			for (uint alt_idx = 0; alt_idx < cur_var_2->numAlts(); alt_idx++) {

				addOverlapStats(&allele_overlap_stats, chromosome_type, nullptr, 0, cur_var_2, alt_idx);
			}

			num_variants_2++;
			delete cur_var_2;
			
			cur_var_2_parsed = vcf_reader_2.getNextVariant(&cur_var_2);
		}

		cout << "[" << Utils::getLocalTime() << "] Finished chromosome " << contig.id() << endl;
	}

	assert(!cur_var_1_parsed);
	assert(!(vcf_reader_1.getNextVariant(&cur_var_1)));

	assert(!cur_var_2_parsed);
	assert(!(vcf_reader_2.getNextVariant(&cur_var_2)));


    ofstream allele_overlap_stats_writer(string(argv[3]) + "_allele_overlap_stats.txt");
    assert(allele_overlap_stats_writer.is_open());

	allele_overlap_stats_writer << "Count\tChromosomeType\tFilters1\tCallProbablity1\tFilters2\tCallProbablity2\tAlleleType\tAlleleLength\tNumAmbiguous\tAlleleSVLength\n";

	for (auto &stats: allele_overlap_stats) {

		allele_overlap_stats_writer << stats.second << "\t" << stats.first << "\n";
	}

	allele_overlap_stats_writer.close();


    ofstream concordance_trio_stats_writer(string(argv[3]) + "_concordance_trio_stats.txt");
    assert(concordance_trio_stats_writer.is_open());

	concordance_trio_stats_writer << "Count\tChromosomeType\tTrioId\tIsConcordant1\tIsConcordant2\tFilter1\tFilter2\tIsReferenceCall1\tIsReferenceCall2\tIsParentsBiAllelicHeterozygote1\tIsParentsBiAllelicHeterozygote2\tMinGPP1\tMinGPP2\tVariantType\tHasMissing\tHasRepeat\tNumAlleles\tMaxAlleleLength\tMaxNumAmbiguous\tMaxAbsAlleleSVLength\n";

	for (auto &stats: concordance_trio_stats) {

		concordance_trio_stats_writer << stats.second << "\t" << stats.first << "\n";
	}

	concordance_trio_stats_writer.close();


	cout << "\n[" << Utils::getLocalTime() << "] Parsed " << num_variants_1 << " variants in callset1 and " <<  num_variants_2 << " variants in callset2" << endl;
    
	cout << "\n[" << Utils::getLocalTime() << "] Found " << num_variants_overlap << " variants identical between callsets" << endl;

    cout << endl;

	return 0;
}
