
/*
Auxiliaries.cpp - This file is part of BayesTyper (v0.9)


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

namespace Auxiliaries {

    uint leftTrimAllelePair(Allele * allele_1, Allele * allele_2, const bool full_trim) {

		assert(!(allele_1->seq().empty()));
		assert(!(allele_2->seq().empty()));

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

        assert(!(allele_1->seq().empty()));
        assert(!(allele_2->seq().empty()));

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

		assert(!(allele_1->seq().empty()));
		assert(!(allele_2->seq().empty()));

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

        assert(!(allele_1->seq().empty()));
        assert(!(allele_2->seq().empty()));

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

        assert(!(main_allele.seq().empty()));
        assert(!(reference_allele.seq().empty()));

        assert(!(reference_allele.isMissing()));

    	if (main_allele.isMissing()) {

    		return AlleleAttributes(Type::Missing, 0, 0, 0);
    	}

    	if (main_allele == reference_allele) {

    		return AlleleAttributes(Type::Reference, main_allele.seq().size(), count(main_allele.seq().begin(), main_allele.seq().end(), 'N'), 0);
    	}

    	Allele trimmed_main_allele = main_allele;
    	Allele trimmed_reference_allele = reference_allele;

    	fullTrimAllelePair(&trimmed_main_allele, &trimmed_reference_allele);
        assert(!(trimmed_main_allele.seq().empty()) or !(trimmed_reference_allele.seq().empty()));

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

        assert(variant.numAlts() > hasMissing(variant));

        if (variant.numAlts() > (1 + static_cast<uint>(hasMissing(variant)))) {

            return "Multi";
        } 

        auto cur_type = alleleAttributes(variant.alt(0), variant.ref());

        assert(cur_type.type != Type::Reference);
        assert(cur_type.type != Type::Missing);

        return cur_type.typeStr();
    }

    bool hasMissing(Variant & variant) {

        assert(!(variant.allele(0).isMissing()));

        for (uint alt_idx = 0; alt_idx < variant.numAlts(); alt_idx++) {

            if (variant.alt(alt_idx).isMissing()) {

                return true;
            }
        }

        return false;
    }

    bool hasRepeat(Variant & variant) {

        for (uint alt_idx = 0; alt_idx < variant.numAlts(); alt_idx++) {

            auto rma_value = variant.alt(alt_idx).info().getValue<string>("RMA");

            if (rma_value.second) {

                if (rma_value.first != ".") {

                    return true;
                }
            }
        }

        return false;
    }

	uint rightTrimVariant(Variant * variant) {

	    uint trimmed = 0;
	    bool has_equal = true;

	    while (has_equal) {

	        if (variant->ref().seq().size() == 1) {

	            has_equal = false;
	            break;
	        }

	        for (uint alt_idx = 0; alt_idx < variant->numAlts(); alt_idx++) {

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

        if (allele.seq().find_first_of("N") != string::npos) {

            return true;
        
        } else {

            return false;
        }
    }
    
    bool hasAnnotatedAlternativeAllele(Variant & variant) {

        for (uint alt_idx = 0; alt_idx < variant.numAlts(); alt_idx++) {

            if (isAlleleAnnotated(variant.alt(alt_idx))) {

                return true;
            }
        }

        return false;
    }

    bool isAlleleAnnotated(Allele & allele) {

        auto annotation = allele.info().getValue<string>("AAI");
        
        if (annotation.second) {

            assert(!(annotation.first.empty()));

            if (annotation.first != ".") {

                return true;
            
            } else {

                return false;
            }

        } else {

            return false;
        }
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

    pair<float, bool> getGenotypePosterior(Sample & sample) {

        if (sample.callStatus() != Sample::CallStatus::Complete) {

            return make_pair(0, false);            
        }

        if (sample.genotypeEstimate().size() == 2) {

            return sample.genotypeInfo().at(sample.twoToOneDimIdx(sample.genotypeEstimate())).getValue<float>("GPP");

        } else if (sample.genotypeEstimate().size() == 1) {

            return sample.genotypeInfo().at(sample.genotypeEstimate().front()).getValue<float>("GPP");
        
        } else {

            return make_pair(0, false);
        }

        return make_pair(0, false);            
    }

    pair<float, bool> getMaxGenotypePosterior(Sample & sample) {

        if (sample.genotypeInfo().empty()) {

            assert(sample.ploidy() == Sample::Ploidy::Zeroploid);
            return make_pair(0, false);            
        }

        assert(sample.ploidy() != Sample::Ploidy::Zeroploid);
        assert(sample.ploidy() != Sample::Ploidy::Polyploid);

        float max_gpp = 0;
        bool has_gpp = false;

        for (auto & info: sample.genotypeInfo()) {

            auto gpp_value = info.getValue<float>("GPP");

            if (gpp_value.second) {

                max_gpp = max(max_gpp, gpp_value.first);
                has_gpp = true;
            }
        }

        assert((max_gpp > 0) or Utils::floatCompare(max_gpp, 0));
        assert((max_gpp < 1) or Utils::floatCompare(max_gpp, 1));
        
        return make_pair(max_gpp, has_gpp);            
    }

    string AlleleAttributes::typeStr() const {

        stringstream type_ss;
        type_ss << type;
        return type_ss.str();
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

    vector<uint> getCalledAlleleIdxsSorted(Sample & sample, const float min_map) {

        vector<uint> called_allele_idxs;
        
        for (uint allele_idx = 0; allele_idx < sample.alleleInfo().size(); allele_idx++) {

            auto map_value = sample.alleleInfo().at(allele_idx).getValue<float>("MAP");
            
            if (map_value.second) {

                if (map_value.first >= min_map) {

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

                if (!(Utils::floatCompare(acp_value.first, 0))) {

                    non_zero_prob_allele_idxs.push_back(allele_idx);
                }
            }
        }

        return non_zero_prob_allele_idxs;
    }

    vector<uint> getNonZeroProbAlleleIdxsSorted(Sample & sample) {

        vector<uint> non_zero_prob_allele_idxs;

        for (uint allele_idx = 0; allele_idx < sample.alleleInfo().size(); allele_idx++) {

            auto map_value = sample.alleleInfo().at(allele_idx).getValue<float>("MAP");

            if (map_value.second) {

                if (!(Utils::floatCompare(map_value.first, 0))) {

                    non_zero_prob_allele_idxs.push_back(allele_idx);
                }
            }
        }

        return non_zero_prob_allele_idxs;
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
