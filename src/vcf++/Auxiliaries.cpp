
/*
Auxiliaries.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#include <iostream>
#include <algorithm>

#include "Auxiliaries.hpp"
#include "Utils.hpp"
#include "Allele.hpp"
#include "Variant.hpp"
#include "JoiningString.hpp"
#include "Stats.hpp"
#include "Contig.hpp"

namespace Auxiliaries {

    uint leftTrimAllelePair(Allele * allele_1, Allele * allele_2, const bool full_trim) {

        assert(!allele_1->isID());
        assert(!allele_2->isID());

		assert(!allele_1->seq().empty());
		assert(!allele_2->seq().empty());

		if (allele_1->isMissing() or allele_2->isMissing()) {

			return 0;
		}

		uint trimmed = 0;

        while ((allele_1->seq().size() > 1) and (allele_2->seq().size() > 1)) {

			if (allele_1->seq().front() != allele_2->seq().front()) {

				break;
			}

			allele_1->seq().erase(allele_1->seq().begin());
			allele_2->seq().erase(allele_2->seq().begin());

			trimmed++;
		}

        assert(!allele_1->seq().empty());
        assert(!allele_2->seq().empty());

        if (full_trim) {

            if (allele_1->seq().front() == allele_2->seq().front()) {

                allele_1->seq().erase(allele_1->seq().begin());
                allele_2->seq().erase(allele_2->seq().begin());

                trimmed++;            
            }
        } 

		return trimmed;
    }

    uint rightTrimAllelePair(Allele * allele_1, Allele * allele_2) {

        assert(!allele_1->isID());
        assert(!allele_2->isID());

		assert(!allele_1->seq().empty());
		assert(!allele_2->seq().empty());

		if (allele_1->isMissing() or allele_2->isMissing()) {

			return 0;
		}

		uint trimmed = 0;

		while ((allele_1->seq().size() > 1) and (allele_2->seq().size() > 1)) {

			if (allele_1->seq().back() != allele_2->seq().back()) {

				break;
			}

			allele_1->seq().pop_back();
			allele_2->seq().pop_back();

			trimmed++;
		}

        assert(!allele_1->seq().empty());
        assert(!allele_2->seq().empty());

		return trimmed;
    }

    pair<uint, uint> partialTrimAllelePair(Allele * allele_1, Allele * allele_2) {

    	pair<uint, uint> trimmed;

    	trimmed.second = rightTrimAllelePair(allele_1, allele_2);
    	trimmed.first = leftTrimAllelePair(allele_1, allele_2, false);

    	return trimmed;
    }

    pair<uint, uint> fullTrimAllelePair(Allele * allele_1, Allele * allele_2) {

        pair<uint, uint> trimmed;

        trimmed.second = rightTrimAllelePair(allele_1, allele_2);
        trimmed.first = leftTrimAllelePair(allele_1, allele_2, true);

        return trimmed;
    }

    AlleleAttributes alleleAttributes(Allele & main_allele, Allele & reference_allele) {

        assert(!main_allele.isID());
        assert(!reference_allele.isID());

        assert(!main_allele.seq().empty());
        assert(!reference_allele.seq().empty());

        assert(!reference_allele.isMissing());

    	if (main_allele.isMissing()) {

    		return AlleleAttributes(Type::Missing, 0, 0, 0);
    	}

    	if (main_allele == reference_allele) {

    		return AlleleAttributes(Type::Reference, main_allele.seq().size(), count(main_allele.seq().begin(), main_allele.seq().end(), 'N'), 0);
    	}

    	Allele trimmed_main_allele = main_allele;
    	Allele trimmed_reference_allele = reference_allele;

    	fullTrimAllelePair(&trimmed_main_allele, &trimmed_reference_allele);
        assert(!trimmed_main_allele.seq().empty() or !trimmed_reference_allele.seq().empty());

        uint trimmed_main_allele_length = trimmed_main_allele.seq().size();
        uint trimmed_reference_allele_length = trimmed_reference_allele.seq().size();

        uint trimmed_main_allele_num_ambiguous = count(trimmed_main_allele.seq().begin(), trimmed_main_allele.seq().end(), 'N');

    	if (trimmed_main_allele_length == trimmed_reference_allele_length) {

            auto allele_type = Type::Complex;

    		if (trimmed_main_allele_length == 1) {

	    		allele_type = Type::SNP;

    		} else if (isInversion(trimmed_main_allele, trimmed_reference_allele, 0.95, 10)) {

                allele_type = Type::Inversion;
    		} 

	    	return AlleleAttributes(allele_type, trimmed_main_allele_length, trimmed_main_allele_num_ambiguous, 0);

    	} else {

            auto allele_type = Type::Complex;

            if (trimmed_main_allele_length == 0) {

                allele_type = Type::Deletion;

            } else if (trimmed_reference_allele_length == 0) {

                allele_type = Type::Insertion;
            } 

            return AlleleAttributes(allele_type, trimmed_main_allele_length, trimmed_main_allele_num_ambiguous, trimmed_main_allele_length - trimmed_reference_allele_length);          
        }
    }

    bool isInversion(Allele & main_allele, Allele & reference_allele, const float min_match_fraction, const uint min_size) {

        assert(!main_allele.isID());
        assert(!reference_allele.isID());

        assert(!reference_allele.isMissing());

        if (main_allele.isMissing()) {

            return false;
        }

    	if (main_allele.seq().size() != reference_allele.seq().size()) {

    		return false;
    	}

        if (main_allele.seq().size() < min_size) {

            return false;
        }

        string main_allele_rv = reverseComplementSequence(main_allele.seq());
        assert(main_allele_rv.size() == reference_allele.seq().size());

    	auto main_rv_it = main_allele_rv.begin();
    	auto reference_rit = reference_allele.seq().begin();

    	uint num_correct_bases = 0;

    	while (main_rv_it != main_allele_rv.end()) {

            if ((*main_rv_it == *reference_rit) and (*main_rv_it != 'N')) {

                num_correct_bases++;                
            }

    		main_rv_it++;
    		reference_rit++;
    	}

    	assert(num_correct_bases <= main_allele_rv.size());
    	assert(reference_rit == reference_allele.seq().end());

    	if ((static_cast<float>(num_correct_bases)/main_allele_rv.size()) < min_match_fraction) {

    		return false;

    	} else {

    		return true;
    	}
    }

    string reverseComplementSequence(const string & sequence) {

        string rv_sequence;
        rv_sequence.reserve(sequence.size());

        for (auto sequence_rit = sequence.rbegin(); sequence_rit != sequence.rend(); sequence_rit++) {

            if ((*sequence_rit == 'A') or (*sequence_rit == 'a')) {

                rv_sequence += "T";

            } else if ((*sequence_rit == 'C') or (*sequence_rit == 'c')) {

                rv_sequence += "G";

            } else if ((*sequence_rit == 'G') or (*sequence_rit == 'g')) {

                rv_sequence += "C";

            } else if ((*sequence_rit == 'T') or (*sequence_rit == 't')) {

                rv_sequence += "A";

            } else {

                assert((*sequence_rit == 'N') or (*sequence_rit == 'n'));
                rv_sequence += "N";
            }
        }   

        assert(rv_sequence.size() == sequence.size());

        return rv_sequence;
    }

    string variantType(Variant & variant) {

        assert(variant.numAlts() > static_cast<uint>(hasMissing(variant)));

        if (variant.numAlts() > (1 + static_cast<uint>(hasMissing(variant)))) {

            return "Multi";
        } 

        auto cur_type = alleleAttributes(variant.alt(0), variant.ref());

        assert(cur_type.type != Type::Reference);
        assert(cur_type.type != Type::Missing);

        return cur_type.typeStr();
    }

    pair<string, bool> variantOrigins(Variant & variant) {

        vector<string> variant_origins;
        
        for (uint alt_idx = 0; alt_idx < variant.numAlts(); alt_idx++) {

            auto aco_value = variant.alt(alt_idx).info().getValue<string>("ACO");

            if (aco_value.second) {

                if (aco_value.first != ".") {

                    assert(!variant.alt(alt_idx).isMissing());

                    auto aco_value_split = Utils::splitString(aco_value.first, ':');
                    assert(find(aco_value_split.begin(), aco_value_split.end(), ".") == aco_value_split.end());

                    for (auto & origin: aco_value_split) {

                        if (find(variant_origins.begin(), variant_origins.end(), origin) == variant_origins.end()) {

                            variant_origins.push_back(origin);
                        }
                    }                           
                }
            }
        }

        if (variant_origins.empty()) {

            return make_pair("", false);
       
        } else {

            sort(variant_origins.begin(), variant_origins.end());

            JoiningString variant_origins_str(':');
            variant_origins_str.join(variant_origins);

            return make_pair(variant_origins_str.str(), true);
        }
    }

    bool hasMissing(Variant & variant) {

        assert(!variant.allele(0).isMissing());

        for (uint alt_idx = 0; alt_idx < variant.numAlts(); alt_idx++) {

            if (variant.alt(alt_idx).isMissing()) {

                return true;
            }
        }

        return false;
    }

    bool isAnnotated(Variant & variant) {

        for (uint alt_idx = 0; alt_idx < variant.numAlts(); alt_idx++) {

            if (isAnnotated(variant.alt(alt_idx))) {

                return true;
            }
        }

        return false;
    }

    bool isAnnotated(Allele & allele) {

        auto aat_value = allele.info().getValue<string>("AAI");
        
        if (aat_value.second) {

            assert(!aat_value.first.empty());

            if (aat_value.first != ".") {

                return true; 
            } 
        } 

        return false;
    }

    bool hasRepeat(Variant & variant) {

        for (uint alt_idx = 0; alt_idx < variant.numAlts(); alt_idx++) {

            if (hasRepeat(variant.alt(alt_idx))) {

                return true;
            }
        }

        return false;
    }

    bool hasRepeat(Allele & allele) {

        auto rma_value = allele.info().getValue<string>("RMA");

        if (rma_value.second) {

            assert(!rma_value.first.empty());

            if (rma_value.first != ".") {

                return true;
            }
        }

        return false;
    }

    uint repeatLength(Allele & allele, const string & name) {

        auto rma_value = allele.info().getValue<string>("RMA");

        if (rma_value.second) {

            if (rma_value.first != ".") {

                auto rma_value_split = Utils::splitString(rma_value.first, ':');

                for (auto & repeat: rma_value_split) {

                    auto repeat_split = Utils::splitString(repeat, '#');
                    assert(repeat_split.size() == 2);

                    if (repeat_split.front() == name) {

                        return stoi(repeat_split.back());
                    }
                }
            }
        }

        return 0;
    }

	uint rightTrimVariant(Variant * variant) {

	    uint trimmed = 0;
	    bool has_equal = true;

	    while (has_equal) {

	        if (variant->ref().seq().size() == 1) {

	            has_equal = false;
	            break;
	        }

            assert(!variant->ref().isID());

	        for (uint alt_idx = 0; alt_idx < variant->numAlts(); alt_idx++) {

                assert(!variant->alt(alt_idx).isID());

	            if (variant->alt(alt_idx).seq().size() == 1) {

	                has_equal = false;
	                break;
	            }

	            if (variant->ref().seq().back() != variant->alt(alt_idx).seq().back()) {

	                has_equal = false;
	                break;
	            }
	        }

	        if (has_equal) {

	            trimmed++;

	            variant->ref().seq().pop_back();
	            assert(variant->ref().seq().size() > 0);

	            for (uint alt_idx = 0; alt_idx < variant->numAlts(); alt_idx++) {

	                variant->alt(alt_idx).seq().pop_back();
	                assert(variant->alt(alt_idx).seq().size() > 0);
	            }
	        }
	    }

	    return trimmed;
	}

    bool hasAmbiguous(Variant & variant) {

        for (uint allele_idx = 0; allele_idx < variant.numAlls(); allele_idx++) {

            if (hasAmbiguous(variant.allele(allele_idx))) {

                return true;
            }
        }

        return false;
    }

    bool hasAmbiguous(Allele & allele) {

        if (allele.isID() or allele.isMissing()) {

            return false;

        } else if (allele.seq().find_first_not_of("ACGT") != string::npos) {

            return true;
        } 

        return false;
    }

    bool isAlleleCalled(Allele & allele, const float min_acp) {

        auto acp = allele.info().getValue<float>("ACP");
        
        if (acp.second) {

            if (acp.first >= min_acp) {

                return true;
            
            } else {

                return false;
            }

        } else {

            return false;
        }
    }    

    pair<float, bool> getMaxGenotypePosterior(Sample & sample) {

        pair<float, bool> max_gpp(0, false);

        auto max_gpp_idx = getMaxGenotypePosteriorIndex(sample);

        if (max_gpp_idx.second) {

            assert(max_gpp_idx.first >= 0);

            max_gpp.first = sample.genotypeInfo().at(max_gpp_idx.first).getValue<float>("GPP").first;
            max_gpp.second = max_gpp_idx.second;

            assert((max_gpp.first > 0) or Utils::floatCompare(max_gpp.first, 0));
            assert((max_gpp.first < 1) or Utils::floatCompare(max_gpp.first, 1));
        }

        return max_gpp;            
    }

    pair<int, bool> getMaxGenotypePosteriorIndex(Sample & sample) {

        pair<int, bool> max_gpp_idx(-1, false);
        float max_gpp_value = 0;

        for (uint gpp_idx = 0; gpp_idx < sample.genotypeInfo().size(); gpp_idx++) {

            auto gpp_value = sample.genotypeInfo().at(gpp_idx).getValue<float>("GPP");

            if (gpp_value.second) {

                if (Utils::floatCompare(max_gpp_value, gpp_value.first)) {

                    max_gpp_idx.first = gpp_idx;
                    max_gpp_idx.second = false;

                } else if (max_gpp_value < gpp_value.first) {

                    max_gpp_idx.first = gpp_idx;
                    max_gpp_idx.second = true;

                    max_gpp_value = gpp_value.first;
                }
            
            } else {

                return make_pair(-1, false);
            }
        }

        return max_gpp_idx;
    }

    void resetFilters(Variant * variant) {

        variant->setFilters({});

        for (auto & sample_id: variant->sampleIds()) {

            Sample * sample = &(variant->getSample(sample_id));

            auto max_gpp_idx = getMaxGenotypePosteriorIndex(*sample);

            if (max_gpp_idx.first >= 0) {

                if (max_gpp_idx.second) {

                    assert(sample->ploidy() != Sample::Ploidy::Zeroploid);
                    assert(sample->ploidy() != Sample::Ploidy::Polyploid);

                    if (sample->ploidy() == Sample::Ploidy::Diploid) {

                        sample->newGenotypeEstimate(sample->oneToTwoDimIdx(max_gpp_idx.first));

                    } else if (sample->ploidy() == Sample::Ploidy::Haploid) {

                        sample->newGenotypeEstimate(vector<ushort>(1, max_gpp_idx.first));
                    }                 
                
                } else {

                    sample->newGenotypeEstimate(vector<ushort>());
                }
            
                for (uint allele_idx = 0; allele_idx < sample->alleleInfo().size(); allele_idx++) {

                    sample->alleleInfo().at(allele_idx).setValue<int>("SAF", 0);
                }
            
            } else {

                sample->newGenotypeEstimate(vector<ushort>());                
            }
        }
    }

    void updateAlleleStatsAndCallProb(Variant * variant) {

        auto allele_stats = Stats::calcAlleleStats(*variant);
        assert(allele_stats.allele_counts.size() == variant->numAlls());
        assert(allele_stats.allele_freqs.size() == variant->numAlls());

        variant->info().setValue<int>("AN", allele_stats.allele_count_sum);

        auto call_probs = Stats::calcCallProbs(*variant);
        assert(call_probs.allele_call_probs.size() == variant->numAlls());

        for (uint allele_idx = 0; allele_idx < variant->numAlls(); allele_idx++) {

            variant->allele(allele_idx).info().setValue<float>("ACP", call_probs.allele_call_probs.at(allele_idx));
    
            if (allele_idx > 0) {

                variant->allele(allele_idx).info().setValue<int>("AC", allele_stats.allele_counts.at(allele_idx));
                variant->allele(allele_idx).info().setValue<float>("AF", allele_stats.allele_freqs.at(allele_idx));
            }
        }

        variant->setQual(make_pair(call_probs.variant_quality, true));
    }

    vector<uint> getCalledAlleleIdxsSorted(Variant & var, const float min_acp) {

        vector<uint> called_allele_idxs;
        
        for (uint allele_idx = 0; allele_idx < var.numAlls(); allele_idx++) {

            if (isAlleleCalled(var.allele(allele_idx), min_acp)) {

                called_allele_idxs.push_back(allele_idx);
            }
        }

        return called_allele_idxs;
    }

    vector<uint> getCalledAlleleIdxsSorted(Sample & sample, const float min_app) {

        vector<uint> called_allele_idxs;
        
        for (uint allele_idx = 0; allele_idx < sample.alleleInfo().size(); allele_idx++) {

            auto app_value = sample.alleleInfo().at(allele_idx).getValue<float>("APP");
            
            if (app_value.second) {

                if (app_value.first >= min_app) {

                    called_allele_idxs.push_back(allele_idx);            
                } 
            }
        }

        return called_allele_idxs;
    }

    vector<uint> getNonZeroProbAlleleIdxsSorted(Variant & var) {

        vector<uint> non_zero_prob_allele_idxs;

        for (uint allele_idx = 0; allele_idx < var.numAlls(); allele_idx++) {

            auto acp_value = var.allele(allele_idx).info().getValue<float>("ACP");

            if (acp_value.second) {

                if (!Utils::floatCompare(acp_value.first, 0)) {

                    non_zero_prob_allele_idxs.push_back(allele_idx);
                }
            }
        }

        return non_zero_prob_allele_idxs;
    }

    vector<uint> getNonZeroProbAlleleIdxsSorted(Sample & sample) {

        vector<uint> non_zero_prob_allele_idxs;

        for (uint allele_idx = 0; allele_idx < sample.alleleInfo().size(); allele_idx++) {

            auto app_value = sample.alleleInfo().at(allele_idx).getValue<float>("APP");

            if (app_value.second) {

                if (!Utils::floatCompare(app_value.first, 0)) {

                    non_zero_prob_allele_idxs.push_back(allele_idx);
                }
            }
        }

        return non_zero_prob_allele_idxs;
    }

    void removeNonRelevantFilterDescriptors(VcfMetaData * meta_data, const unordered_set<string> & relevant_filter_discriptors) {

        auto filter_it = meta_data->filterDescriptors().begin();

        while (filter_it != meta_data->filterDescriptors().end()) {

            if (relevant_filter_discriptors.count(filter_it->first) > 0) {

                filter_it++;

            } else {

                filter_it = meta_data->filterDescriptors().erase(filter_it);
            }
        }
    }

    void removeNonRelevantInfoDescriptors(VcfMetaData * meta_data, const unordered_set<string> & relevant_info_discriptors) {

        auto info_it = meta_data->infoDescriptors().begin();

        while (info_it != meta_data->infoDescriptors().end()) {

            if (relevant_info_discriptors.count(info_it->first) > 0) {

                info_it++;

            } else {

                info_it = meta_data->infoDescriptors().erase(info_it);
            }
        }
    }

    void removeNonRelevantFormatDescriptors(VcfMetaData * meta_data, const unordered_set<string> & relevant_format_discriptors) {

        auto format_it = meta_data->formatDescriptors().begin();

        while (format_it != meta_data->formatDescriptors().end()) {

            if (relevant_format_discriptors.count(format_it->first) > 0) {

                format_it++;

            } else {

                format_it = meta_data->formatDescriptors().erase(format_it);
            }
        }
    }

    string AlleleAttributes::typeStr() const {

        stringstream type_ss;
        type_ss << type;
        return type_ss.str();
    }

    pair<uint, string> getHomopolymerInfo(uint start_pos, const string & seq) { // start_pos zero_based

    	const char start_nt = seq.at(start_pos);

    	uint hpl_start_idx = start_pos;
 
    	while ((0 <= hpl_start_idx) and (seq.at(hpl_start_idx) == start_nt)) {

    		--hpl_start_idx;
    	}

    	++hpl_start_idx;

    	uint hpl_end_idx = start_pos;

    	while ((hpl_end_idx < seq.size()) and (seq.at(hpl_end_idx) == start_nt)) {

    		++hpl_end_idx;
    	}

    	--hpl_end_idx;

    	assert(0 <= hpl_start_idx);
        assert(hpl_end_idx < seq.size());

        assert(hpl_start_idx <= hpl_end_idx);

    	return make_pair(hpl_end_idx - hpl_start_idx + 1, seq.substr(start_pos, 1));
    }

    vector<bool> getHomopolymerAlleles(Variant & cur_var, const string & homopolymer_base, const uint max_length_difference) {

        vector<bool> homopolymer_alleles(cur_var.numAlls(), false);

        for (uint allele_idx = 0; allele_idx < cur_var.numAlls(); allele_idx++) {

            for (uint homopolymer_allele_idx = allele_idx + 1; homopolymer_allele_idx < cur_var.numAlls(); homopolymer_allele_idx++) {

                if (homopolymer_alleles.at(allele_idx) and homopolymer_alleles.at(homopolymer_allele_idx)) {

                    continue;
                }

                Allele cur_ref_allele_1 = cur_var.ref();
                Allele cur_ref_allele_2 = cur_var.ref();

                Allele cur_allele_1 = cur_var.allele(allele_idx);
                Allele cur_allele_2 = cur_var.allele(homopolymer_allele_idx);

                if (cur_allele_1 != cur_allele_2) {

                    Auxiliaries::fullTrimAllelePair(&cur_ref_allele_1, &cur_allele_1);
                    Auxiliaries::fullTrimAllelePair(&cur_ref_allele_2, &cur_allele_2);

                    bool is_cur_ref_allele_1_homopolymer = cur_ref_allele_1.seq().find_first_not_of(homopolymer_base) == string::npos;
                    bool is_cur_ref_allele_2_homopolymer = cur_ref_allele_2.seq().find_first_not_of(homopolymer_base) == string::npos;
                    
                    bool is_cur_allele_1_homopolymer = cur_allele_1.seq().find_first_not_of(homopolymer_base) == string::npos;
                    bool is_cur_allele_2_homopolymer = cur_allele_2.seq().find_first_not_of(homopolymer_base) == string::npos;

                    if (is_cur_ref_allele_1_homopolymer and is_cur_ref_allele_2_homopolymer and is_cur_allele_1_homopolymer and is_cur_allele_2_homopolymer) {

                        uint cur_ref_allele_length_diff = abs(static_cast<int>(cur_ref_allele_1.seq().size()) - static_cast<int>(cur_ref_allele_2.seq().size()));
                        uint cur_allele_length_diff = abs(static_cast<int>(cur_allele_1.seq().size()) - static_cast<int>(cur_allele_2.seq().size()));

                        assert((cur_ref_allele_length_diff + cur_allele_length_diff) > 0);

                        if ((cur_ref_allele_length_diff + cur_allele_length_diff) <= max_length_difference) {

                            assert(!cur_var.allele(allele_idx).isMissing());
                            homopolymer_alleles.at(allele_idx) = true;     

                            assert(!cur_var.allele(homopolymer_allele_idx).isMissing());
                            homopolymer_alleles.at(homopolymer_allele_idx) = true;              
                        }
                    } 
                }
            }
        }

        return homopolymer_alleles;
    }

    vector<Contig> mergeContigs(const vector<Contig> & contigs_1, const vector<Contig> & contigs_2) {

        bool incorrect_order = false;

        vector<Contig> contigs_merged;
        contigs_merged.reserve(contigs_1.size() + contigs_2.size());

        auto contigs_1_it = contigs_1.begin();
        auto contigs_2_bit = contigs_2.begin();

        while (contigs_1_it != contigs_1.end()) {

            if (incorrect_order or find(contigs_merged.begin(), contigs_merged.end(), *contigs_1_it) != contigs_merged.end()) {

                incorrect_order = true;
                break;
            }

            auto contigs_2_eit = find(contigs_2.begin(), contigs_2.end(), *contigs_1_it);

            if (contigs_2_eit != contigs_2.end()) {

                while (contigs_2_bit != contigs_2_eit) {

                    if (find(contigs_merged.begin(), contigs_merged.end(), *contigs_2_bit) != contigs_merged.end()) {

                        incorrect_order = true;
                        break;
                    }

                    contigs_merged.push_back(*contigs_2_bit);
                    contigs_2_bit++;
                }

                contigs_2_bit++;
            }

            contigs_merged.push_back(*contigs_1_it);
            contigs_1_it++;
        }

        while (contigs_2_bit != contigs_2.end()) {

            if (find(contigs_merged.begin(), contigs_merged.end(), *contigs_2_bit) != contigs_merged.end()) {

                incorrect_order = true;
                break;
            }

            contigs_merged.push_back(*contigs_2_bit);
            contigs_2_bit++;
        }

        if (incorrect_order) {

            cerr << "\nERROR: Contigs should be sorted in the same order across VCF files (both header and variants)\n" << endl;
            exit(1);
        }

        assert(contigs_merged.size() <= (contigs_1.size() + contigs_2.size()));
        return contigs_merged;
    }
}


ostream& operator<< (ostream & os, Auxiliaries::Type type) {

    switch (type) {

        case Auxiliaries::Type::Reference:

            return os << "Reference";

        case Auxiliaries::Type::SNP:

            return os << "SNP";

        case Auxiliaries::Type::Insertion:

            return os << "Insertion";

        case Auxiliaries::Type::Deletion:

            return os << "Deletion";

        case Auxiliaries::Type::Inversion:

            return os << "Inversion";

        case Auxiliaries::Type::Complex:

            return os << "Complex";

        case Auxiliaries::Type::Missing:

            return os << "Missing";

        default:

        	assert(false);
    }

    return os;
}

bool operator==(const Auxiliaries::AlleleAttributes & lhs, const Auxiliaries::AlleleAttributes & rhs) {
    
    return ((lhs.type == rhs.type) and (lhs.length == rhs.length) and (lhs.num_ambiguous == rhs.num_ambiguous) and (lhs.sv_length == rhs.sv_length));
}

bool operator!=(const Auxiliaries::AlleleAttributes & lhs, const Auxiliaries::AlleleAttributes & rhs) {
    
    return !(lhs == rhs);
}
