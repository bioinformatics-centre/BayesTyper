
/*
getSummary.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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
#include <iostream>
#include <math.h>
#include <unordered_set>
#include <algorithm>

#include "VcfFile.hpp"
#include "Variant.hpp"
#include "JoiningString.hpp"
#include "Auxiliaries.hpp"
#include "Stats.hpp"
#include "Trio.hpp"


static const vector<string> variant_attributes({"Count","ChromType","Filter","VariantType","HasMissing","HasRedundant","NumAlleles","EffectiveNumAlleles","MaxAltACP","MaxAltAC","AN","ACO","HPL","HasHomopolymer","HTV","NumCompleteSamples","NumCONCTrue","NumCONCFalse","BASE","CALL","GTCO","MED"});

static const vector<string> allele_attributes({"Count","ChromType","Filter","AlleleType","AlleleLength","AlleleSVLength","IsRedundant","NumAlleles","EffectiveNumAlleles","ACP","AC","AN","ACO","HPL","IsHomopolymer","HTV","NumCompleteSamples","NumCONCTrue","NumCONCFalse","BASE","CALL","GTCO","MED","MinNAK","MinFAK"});


pair<float, bool> parseSampleAlleleValue(const pair<float, bool> & att_value) {

	if (att_value.second) {

		if (Utils::floatCompare(att_value.first, -1)) {

			return make_pair(0, false);

		} else {

			assert(att_value.first >= 0);
			return att_value;
		}

	} else {

		return make_pair(0, false);
	}
}

pair<string, string> parseValuePair(const pair<string, bool> & att_value) {

	if (att_value.second) {

		auto att_value_split = Utils::splitString(att_value.first, ':');
		assert(att_value_split.size() == 2);

		return make_pair(att_value_split.front(), att_value_split.back());

	} else {

		return make_pair("NA", "NA");
	}
}

string convertValueToString(const pair<string, bool> & att_value) {

	if (att_value.second) {

		if (att_value.first == ".") {

			return "NA";

		} else {

			if (att_value.first == "TRUE") {

				return "1";
			
			} else if (att_value.first == "FALSE") {

				return "0";

			} else {

				return att_value.first;
			}
		}

	} else {

		return "NA";
	}
}

string convertValueToString(const pair<bool, bool> & att_value) {

	if (att_value.second) {

		return Utils::boolToString(att_value.first);

	} else {

		return "NA";
	}
}

string convertValueToString(const pair<int, bool> & att_value) {

	if (att_value.second) {

		return to_string(att_value.first);

	} else {

		return "NA";
	}
}

string convertValueToString(const pair<uint, bool> & att_value) {

	if (att_value.second) {

		return to_string(att_value.first);

	} else {

		return "NA";
	}
}

string convertValueToString(const pair<float, bool> & att_value, const uint precision) {

	if (att_value.second) {

		return Utils::floatToString(att_value.first, precision);

	} else {

		return "NA";
	}
}

template<typename ValueType>
pair<ValueType, bool> getPairValue(const pair<ValueType, bool> & first_value, const pair<ValueType, bool> & second_value, const string & value_func) {

	if (first_value.second and second_value.second) {

		if (value_func == "min") {

			return make_pair(min(first_value.first, second_value.first), true);

		} else {

			assert(value_func == "max");
			return make_pair(max(first_value.first, second_value.first), true);
		}
	
	} else if (first_value.second) {

		return first_value;

	} else if (second_value.second) {

		return second_value;
	}

	return make_pair(0, false);
}

void writeSummaryStats(const unordered_map<string, uint> & summary_stats, const string & output_name, const string & header) {

    ofstream stats_outfile(output_name);
    assert(stats_outfile.is_open());

    stats_outfile << header << endl;

	for (auto & line: summary_stats) {

		stats_outfile << line.second << "\t" << line.first << endl;
	}

	stats_outfile.close();
}



int main(int argc, char const *argv[]) {

    if (argc != 3) {

        std::cout << "USAGE: getSummary <variant_file> <output_prefix>" << std::endl;
        return 1;
    }

	cout << "[" << Utils::getLocalTime() << "] Running BayesTyperTools (" << BT_VERSION << ") getSummary script ...\n" << endl;

	GenotypedVcfFileReader vcf_reader(argv[1], false);

	Variant * cur_var;

	uint num_variants = 0;
	uint num_alleles = 0;

	unordered_map<string, uint> variant_summary_stats;
	unordered_map<string, uint> allele_summary_stats;

	while (vcf_reader.getNextVariant(&cur_var)) {

		num_variants++;
		num_alleles += cur_var->numAlls();

		string chrom_type = vcf_reader.metaData().getContig(cur_var->chrom()).typeStr();

        JoiningString filter_str(';');

        if (cur_var->filters().empty()) {

        	filter_str.join(".");

        } else {

            filter_str.join(cur_var->filters());
        }

        auto allele_stats = Stats::calcAlleleStats(*cur_var);
        assert(allele_stats.allele_counts.size() == cur_var->numAlls());
        assert(allele_stats.allele_freqs.size() == cur_var->numAlls());

        auto call_probs = Stats::calcCallProbs(*cur_var);
        assert(call_probs.allele_call_probs.size() == cur_var->numAlls());

        uint effective_num_alleles = 0;
        float max_alt_acp = 0;
        uint max_alt_ac = 0;

		for (uint allele_idx = 0; allele_idx < cur_var->numAlls(); allele_idx++) {

            if (allele_stats.allele_counts.at(allele_idx) > 0) {

                effective_num_alleles++;
            }
    
            if (allele_idx > 0) {

            	max_alt_acp = max(max_alt_acp, call_probs.allele_call_probs.at(allele_idx));
            	max_alt_ac = max(max_alt_ac, allele_stats.allele_counts.at(allele_idx));
            }
        }

        pair<float, bool> ibc_value_1(0, false);
        pair<int, bool> ibc_value_2(0, false);

        auto ibc_values = parseValuePair(cur_var->info().getValue<string>("IBC"));

        if (ibc_values.first != "NA") {

        	ibc_value_1.first = stof(ibc_values.first);
        	ibc_value_1.second = true;

            assert(ibc_values.second != "NA");

            ibc_value_2.first = stoi(ibc_values.second);
            ibc_value_2.second = true;
        }

        pair<int, bool> hpl_value_1(0, false);
        pair<vector<bool>, bool> homopolymer_alleles = make_pair(vector<bool>(cur_var->numAlls(), false), false);

        auto hpl_values = parseValuePair(cur_var->info().getValue<string>("HPL"));

        if (hpl_values.first != "NA") {

            hpl_value_1.first = stoi(hpl_values.first);
            hpl_value_1.second = true;

            assert(hpl_values.second != "NA");

	        homopolymer_alleles.first = Auxiliaries::getHomopolymerAlleles(*cur_var, hpl_values.second, 1);
	        homopolymer_alleles.second = true;
        }

        bool has_homopolymer_allele = (find(homopolymer_alleles.first.begin(), homopolymer_alleles.first.end(), true) != homopolymer_alleles.first.end());
        assert(homopolymer_alleles.first.size() == cur_var->numAlls());

        bool has_redundant = false;
		vector<bool> is_allele_redundant(cur_var->numAlls(), false);

	    for (uint first_allele_idx = 0; first_allele_idx < cur_var->numAlls(); first_allele_idx++) {

	    	if (is_allele_redundant.at(first_allele_idx)) {

	    		break;
	    	}

	        for (uint second_allele_idx = first_allele_idx + 1; second_allele_idx < cur_var->numAlls(); second_allele_idx++) {

	            if (cur_var->allele(first_allele_idx) == cur_var->allele(second_allele_idx)) {

	            	has_redundant = true;
	                is_allele_redundant.at(first_allele_idx) = true;
	                is_allele_redundant.at(second_allele_idx) = true;
	            }
	        }
	    }

	    uint num_complete_samples = 0;
        uint num_conc_true = 0;
        uint num_conc_false = 0;

        pair<string, bool> gtco_value("", false);
        pair<float, bool> med_value(0, false);

        vector<pair<float, bool> > allele_min_nak(cur_var->numAlls(), make_pair(0, false));
        vector<pair<float, bool> > allele_min_fak(cur_var->numAlls(), make_pair(0, false));
        
		for (auto & sample_id: vcf_reader.metaData().sampleIds()) {

			Sample & cur_sample = cur_var->getSample(sample_id);
			assert(cur_sample.ploidy() != Sample::Ploidy::Polyploid);

			auto conc_value = convertValueToString(cur_sample.info().getValue<string>("CONC"));

			if (conc_value == "1") {

				num_conc_true++;

			} else if (conc_value == "0") {

				num_conc_false++;

			} else {

				assert(conc_value == "NA");
			}

			if (cur_sample.callStatus() == Sample::CallStatus::Complete) {

				num_complete_samples++;
			}

			gtco_value = cur_sample.info().getValue<string>("GTCO");
			med_value = cur_sample.info().getValue<float>("MED");

			assert(cur_sample.alleleInfo().size() == cur_var->numAlls());

			for (uint allele_idx = 0; allele_idx < cur_var->numAlls(); allele_idx++) {

				allele_min_nak.at(allele_idx) = getPairValue<float>(allele_min_nak.at(allele_idx), parseSampleAlleleValue(cur_sample.alleleInfo().at(allele_idx).getValue<float>("NAK")), "min"); 
				allele_min_fak.at(allele_idx) = getPairValue<float>(allele_min_fak.at(allele_idx), parseSampleAlleleValue(cur_sample.alleleInfo().at(allele_idx).getValue<float>("FAK")), "min"); 
			}
		}

		if (vcf_reader.metaData().sampleIds().size() != 1) {

       		gtco_value = make_pair("", false);
        	med_value = make_pair(0, false);
		}

		for (uint allele_idx = 0; allele_idx < cur_var->numAlls(); allele_idx++) {
    
			auto allele_attributes = Auxiliaries::alleleAttributes(cur_var->allele(allele_idx), cur_var->ref());

			JoiningString allele_summary_stats_str('\t');

			allele_summary_stats_str.join(convertValueToString(make_pair(chrom_type, true)));
			allele_summary_stats_str.join(convertValueToString(make_pair(filter_str.str(), true)));
			allele_summary_stats_str.join(convertValueToString(make_pair(allele_attributes.typeStr(), true)));
			allele_summary_stats_str.join(convertValueToString(make_pair(allele_attributes.length, true)));
			allele_summary_stats_str.join(convertValueToString(make_pair(allele_attributes.sv_length, true)));
			allele_summary_stats_str.join(convertValueToString(make_pair(static_cast<bool>(is_allele_redundant.at(allele_idx)), true)));
			allele_summary_stats_str.join(convertValueToString(make_pair(cur_var->numAlls(), true)));
			allele_summary_stats_str.join(convertValueToString(make_pair(effective_num_alleles, true)));			
			allele_summary_stats_str.join(convertValueToString(make_pair(call_probs.allele_call_probs.at(allele_idx), true), 2));
			allele_summary_stats_str.join(convertValueToString(make_pair(allele_stats.allele_counts.at(allele_idx), true)));
			allele_summary_stats_str.join(convertValueToString(make_pair(allele_stats.allele_count_sum, true)));
			allele_summary_stats_str.join(convertValueToString(cur_var->allele(allele_idx).info().getValue<string>("ACO")));
			allele_summary_stats_str.join(convertValueToString(hpl_value_1));
			allele_summary_stats_str.join(convertValueToString(make_pair(static_cast<bool>(homopolymer_alleles.first.at(allele_idx)), homopolymer_alleles.second)));
			allele_summary_stats_str.join(convertValueToString(cur_var->info().getValue<string>("HTV")));
			allele_summary_stats_str.join(convertValueToString(make_pair(num_complete_samples, true)));
			allele_summary_stats_str.join(convertValueToString(make_pair(num_conc_true, true)));
			allele_summary_stats_str.join(convertValueToString(make_pair(num_conc_false, true)));
			allele_summary_stats_str.join(convertValueToString(cur_var->info().getValue<string>("BASE")));
			allele_summary_stats_str.join(convertValueToString(cur_var->info().getValue<string>("CALL")));
			allele_summary_stats_str.join(convertValueToString(gtco_value));
			allele_summary_stats_str.join(convertValueToString(med_value, 2));
			allele_summary_stats_str.join(convertValueToString(allele_min_nak.at(allele_idx), 1));
			allele_summary_stats_str.join(convertValueToString(allele_min_fak.at(allele_idx), 2));

			auto allele_summary_stats_it = allele_summary_stats.emplace(allele_summary_stats_str.str(), 0);
			allele_summary_stats_it.first->second++;
		}

		JoiningString variant_summary_stats_str('\t');

		variant_summary_stats_str.join(convertValueToString(make_pair(chrom_type, true)));
		variant_summary_stats_str.join(convertValueToString(make_pair(filter_str.str(), true)));
		variant_summary_stats_str.join(convertValueToString(make_pair(Auxiliaries::variantType(*cur_var), true)));
		variant_summary_stats_str.join(convertValueToString(make_pair(Auxiliaries::hasMissing(*cur_var), true)));
		variant_summary_stats_str.join(convertValueToString(make_pair(has_redundant, true)));
		variant_summary_stats_str.join(convertValueToString(make_pair(cur_var->numAlls(), true)));
		variant_summary_stats_str.join(convertValueToString(make_pair(effective_num_alleles, true)));
		variant_summary_stats_str.join(convertValueToString(make_pair(max_alt_acp, true), 2));
		variant_summary_stats_str.join(convertValueToString(make_pair(max_alt_ac, true)));
		variant_summary_stats_str.join(convertValueToString(make_pair(allele_stats.allele_count_sum, true)));
		variant_summary_stats_str.join(convertValueToString(Auxiliaries::variantOrigins(*cur_var)));
		variant_summary_stats_str.join(convertValueToString(hpl_value_1));
		variant_summary_stats_str.join(convertValueToString(make_pair(has_homopolymer_allele, homopolymer_alleles.second)));
		variant_summary_stats_str.join(convertValueToString(cur_var->info().getValue<string>("HTV")));
		variant_summary_stats_str.join(convertValueToString(make_pair(num_complete_samples, true)));
		variant_summary_stats_str.join(convertValueToString(make_pair(num_conc_true, true)));
		variant_summary_stats_str.join(convertValueToString(make_pair(num_conc_false, true)));
		variant_summary_stats_str.join(convertValueToString(cur_var->info().getValue<string>("BASE")));
		variant_summary_stats_str.join(convertValueToString(cur_var->info().getValue<string>("CALL")));
		variant_summary_stats_str.join(convertValueToString(gtco_value));
		variant_summary_stats_str.join(convertValueToString(med_value, 2));

		auto variant_summary_stats_it = variant_summary_stats.emplace(variant_summary_stats_str.str(), 0);
		variant_summary_stats_it.first->second++;

		delete cur_var;

		if ((num_variants % 100000) == 0) {

			std::cout << "[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants" << endl;
		}
	}


	JoiningString variant_summary_stats_header('\t');
	variant_summary_stats_header.join(variant_attributes);

	writeSummaryStats(variant_summary_stats, string(argv[2]) + "_variant.txt", variant_summary_stats_header.str());

	JoiningString allele_summary_stats_header('\t');
	allele_summary_stats_header.join(allele_attributes);

	writeSummaryStats(allele_summary_stats, string(argv[2]) + "_allele.txt", allele_summary_stats_header.str());

	cout << "\n[" << Utils::getLocalTime() << "] Wrote summary statistics for a total of " << num_variants << " variants and " << num_alleles << " alleles" << endl;
	cout << endl;
}

