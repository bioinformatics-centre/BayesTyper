
/*
getSummary.cpp - This file is part of BayesTyper (v0.9)


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

#include "vcf++/JoiningString.hpp"
#include "vcf++/Auxiliaries.hpp"
#include "vcf++/Stats.hpp"

#include "getSummary.hpp"

namespace GetSummary {

	static const vector<string> variant_misc_attributes({"Count","ChromosomeType","Filter","HasAmbiguous","HasAnnotatedAltAllele","VariantType","NumAlleles","NumCalledAlleles","NumFilteredSamples","NumCalledSamples"});

	static const vector<string> allele_misc_attributes({"Count","ChromosomeType","Filter","IsAnnotated","Origin","IsCalled","AlleleType","AlleleLength","NumAmbiguous","AlleleSVLength","RepeatMaskerAnnotation"});

	static const vector<string> mac_misc_attributes({"Count","ChromosomeType","Filter","AltAlleleHasAnnotation","AltAlleleOrigin","AltAlleleRepeatMaskerAnnotation","AltAlleleType","AltAlleleLength","NumAmbiguous","AltAlleleSVLength","NumCalledParents","MinAlleleCount","MinAlleleCountIsAltAllele"});

	static const vector<string> de_novo_misc_attributes({"Count","Sample","AlleleType","AlleleLength","AlleleSVLength","NumAlleles","NumNonZeroProbAlleles","RepeatMaskerAnnotation"});

	string getAlleleStringAttribute(Allele & allele, const string attribute) {

		auto att_value = allele.info().getValue<string>(attribute);
		
		if (att_value.second) {									
		
			if (att_value.first == ".") {

				return "NoValue";
			
			} else {

				assert(!(allele.isMissing()));
				return att_value.first;
			}

		} else {

			return "Reference";
		}
	}

	void writeSummaryStats(const unordered_map<string, uint> & summary_stats, const string & output_name, const string & header) {

	    ofstream summary_stats_writer(output_name);
	    assert(summary_stats_writer.is_open());
	    
	    summary_stats_writer << header << endl;

		for (auto &line: summary_stats) {

			summary_stats_writer << line.second << "\t" << line.first << endl;
		}

		summary_stats_writer.close();
	}

	void getSummary(const string & vcf_filename, const string & output_prefix, const float min_called_probability, const string & parents_trio_regular_expression, const vector<string> & excluded_sample_ids) {

		cout << "[" << Utils::getLocalTime() << "] Running BayesTyperUtils (" << BTU_VERSION << ") getSummary ...\n" << endl;

		GenotypedVcfFileReader vcf_reader(vcf_filename, true);

        for (auto &sample_id: excluded_sample_ids) {

        	vcf_reader.metaData().rmSample(sample_id);
        }

        cout << "[" << Utils::getLocalTime() << "] Excluded " << excluded_sample_ids.size() << " sample(s)" << endl;

		uint num_parent_samples = 0;
		regex parent_sample_id_regex(parents_trio_regular_expression);

    	for (auto &sample_id: vcf_reader.metaData().sampleIds()) {

        	if (regex_match(sample_id, parent_sample_id_regex)) {

        		num_parent_samples++;
        	}
        }

		cout << "[" << Utils::getLocalTime() << "] Calculating population statistics on " << num_parent_samples << " parent sample(s)\n" << endl;

		Variant * cur_var;

		uint num_variants = 0;
		uint num_alleles = 0;

		uint num_genotypes = 0;

		unordered_map<string, uint> variant_summary_stats;
		unordered_map<string, uint> allele_summary_stats;
		unordered_map<string, uint> mac_summary_stats;
		unordered_map<string, uint> de_novo_summary_stats;
	 
		while (vcf_reader.getNextVariant(&cur_var)) {

			num_variants++;
			num_alleles += cur_var->numAlls();

			assert(cur_var->filters().size() == 1);

			string chrom_type = vcf_reader.metaData().getContig(cur_var->chrom()).typeStr();

			uint num_filtered_samples = 0;
			uint num_filtered_parent_samples = 0;

			uint num_called_samples = 0;
			uint num_called_parent_samples = 0;

			for (auto &sample_id: vcf_reader.metaData().sampleIds()) {

    			num_genotypes++;
    			Sample & cur_sample = cur_var->getSample(sample_id);

				assert(cur_sample.ploidy() != Sample::Ploidy::Polyploid);
				assert(cur_sample.callStatus() != Sample::CallStatus::Partial);

    			if (!(cur_sample.isInformative())) {

    				if (cur_sample.callStatus() == Sample::CallStatus::Missing) {

    					num_filtered_samples++;

    					if (regex_match(sample_id, parent_sample_id_regex)) {

    						num_filtered_parent_samples++;
    					}

    				} else {

    					assert(chrom_type != "Autosomal");
						assert(cur_sample.ploidy() == Sample::Ploidy::Zeroploid);		
    				}

    			} else {

    				assert(cur_sample.callStatus() == Sample::CallStatus::Complete);

    				auto genotype_gpp = Auxiliaries::getGenotypePosterior(cur_sample);
            		assert(genotype_gpp == Auxiliaries::getMaxGenotypePosterior(cur_sample));

    				assert(genotype_gpp.second);

					if (genotype_gpp.first >= min_called_probability) {

    					num_called_samples++;

    					if (regex_match(sample_id, parent_sample_id_regex)) {

    						num_called_parent_samples++;
    					}	
					}
    			}
			}

			assert(num_filtered_samples <= vcf_reader.metaData().sampleIds().size());
			assert(num_filtered_parent_samples <= num_parent_samples);

			auto allele_call_prob = Stats::calcAlleleCallProbAndQualFromAllelePosteriors(cur_var);

			for (uint i = 0; i < cur_var->numAlls(); i++) {

            	assert(!cur_var->allele(i).info().setValue<float>("ACP", allele_call_prob.first.at(i)));
		
				JoiningString allele_summary_stats_str('\t');
				allele_summary_stats_str.join(chrom_type);
				allele_summary_stats_str.join(cur_var->filters().front());
				allele_summary_stats_str.join(Utils::boolToString(Auxiliaries::isAlleleAnnotated(cur_var->allele(i))));
				allele_summary_stats_str.join(getAlleleStringAttribute(cur_var->allele(i), "ACO"));
				allele_summary_stats_str.join(Utils::boolToString(Auxiliaries::isAlleleCalled(cur_var->allele(i), min_called_probability)));

				auto allele_attributes = Auxiliaries::alleleAttributes(cur_var->allele(i), cur_var->ref());
				auto at_value = cur_var->allele(i).info().getValue<string>("AT");
				
				if (at_value.second) {

					assert(allele_attributes.typeStr() == at_value.first);

				} else {

					assert(allele_attributes.type == Auxiliaries::Type::Reference);						
				}

				allele_summary_stats_str.join(allele_attributes.typeStr());
				allele_summary_stats_str.join(to_string(allele_attributes.length));
				allele_summary_stats_str.join(to_string(allele_attributes.num_ambiguous));
				allele_summary_stats_str.join(to_string(allele_attributes.sv_length));
				allele_summary_stats_str.join(getAlleleStringAttribute(cur_var->allele(i), "RMA"));

				auto allele_summary_stats_emplace = allele_summary_stats.emplace(allele_summary_stats_str.str(), 0);
				allele_summary_stats_emplace.first->second++;
			}

			JoiningString variant_summary_stats_str('\t');

			variant_summary_stats_str.join(chrom_type);
			variant_summary_stats_str.join(cur_var->filters().front());
			variant_summary_stats_str.join(Utils::boolToString(Auxiliaries::hasAmbiguous(*cur_var)));				
			variant_summary_stats_str.join(Utils::boolToString(Auxiliaries::hasAnnotatedAlternativeAllele(*cur_var)));				
			variant_summary_stats_str.join(Auxiliaries::variantType(*cur_var));
			variant_summary_stats_str.join(to_string(cur_var->numAlls()));
			variant_summary_stats_str.join(to_string(Auxiliaries::getCalledAlleleIdxsSorted(*cur_var, min_called_probability).size()));
			variant_summary_stats_str.join(to_string(num_filtered_samples));
			variant_summary_stats_str.join(to_string(num_called_samples));

			auto variant_summary_stats_emplace = variant_summary_stats.emplace(variant_summary_stats_str.str(), 0);
			variant_summary_stats_emplace.first->second++;

			if ((cur_var->numAlls() == 2) and (num_filtered_parent_samples == 0)) {

				assert(!(cur_var->allele(1).isMissing()));

				JoiningString mac_summary_stats_str('\t');

				mac_summary_stats_str.join(chrom_type);
				mac_summary_stats_str.join(cur_var->filters().front());
				mac_summary_stats_str.join(Utils::boolToString(Auxiliaries::hasAnnotatedAlternativeAllele(*cur_var)));				
				mac_summary_stats_str.join(getAlleleStringAttribute(cur_var->allele(1), "ACO"));
				mac_summary_stats_str.join(getAlleleStringAttribute(cur_var->allele(1), "RMA"));

				auto allele_attributes = Auxiliaries::alleleAttributes(cur_var->allele(1), cur_var->ref());
				auto at_value = cur_var->allele(1).info().getValue<string>("AT");
				
				assert(at_value.second);
				assert(allele_attributes.typeStr() == at_value.first);

				mac_summary_stats_str.join(allele_attributes.typeStr());
				mac_summary_stats_str.join(to_string(allele_attributes.length));
				mac_summary_stats_str.join(to_string(allele_attributes.num_ambiguous));
				mac_summary_stats_str.join(to_string(allele_attributes.sv_length));
				mac_summary_stats_str.join(to_string(num_called_parent_samples));

				auto allele_stats = Stats::calcAlleleStats(cur_var, parent_sample_id_regex);				
				assert(allele_stats.second);
				
				assert(allele_stats.first.allele_count_sum > 0);
				assert(allele_stats.first.allele_counts.size() == 2);

				if (allele_stats.first.allele_counts.front() <= allele_stats.first.allele_counts.back()) {

					mac_summary_stats_str.join(to_string(allele_stats.first.allele_counts.front()));
					mac_summary_stats_str.join(Utils::boolToString(false));
				
				} else {

					mac_summary_stats_str.join(to_string(allele_stats.first.allele_counts.back()));	
					mac_summary_stats_str.join(Utils::boolToString(true));
				}

				auto mac_summary_stats_emplace = mac_summary_stats.emplace(mac_summary_stats_str.str(), 0);
				mac_summary_stats_emplace.first->second++;
			}

			if (cur_var->info().hasValue("DNE")) {

    			assert(chrom_type == "Autosomal");

				auto dne_att_split = Utils::splitString(cur_var->info().getValue<string>("DNE").first, ',');
				assert(!(dne_att_split.empty()));

				for (auto &dne: dne_att_split) {

					auto dne_value_split = Utils::splitString(dne, ':');
					assert(dne_value_split.size() == 5);

					JoiningString de_novo_summary_stats_str('\t');
					de_novo_summary_stats_str.join(dne_value_split.front());

					assert(stoi(dne_value_split.at(2)) > 0);
					assert(stoi(dne_value_split.at(2)) != stoi(dne_value_split.at(1)));

					auto allele_attributes = Auxiliaries::alleleAttributes(cur_var->allele(stoi(dne_value_split.at(2))), cur_var->allele(stoi(dne_value_split.at(1))));
					assert(allele_attributes.typeStr() == dne_value_split.at(3));
					assert(allele_attributes.num_ambiguous == 0);

					de_novo_summary_stats_str.join(allele_attributes.typeStr());
					de_novo_summary_stats_str.join(to_string(allele_attributes.length));
					de_novo_summary_stats_str.join(to_string(allele_attributes.sv_length));

					de_novo_summary_stats_str.join(to_string(cur_var->numAlls()));
					de_novo_summary_stats_str.join(to_string(Auxiliaries::getNonZeroProbAlleleIdxsSorted(*cur_var).size()));

					if (stoi(dne_value_split.at(1)) == 0) {

						de_novo_summary_stats_str.join(getAlleleStringAttribute(cur_var->allele(stoi(dne_value_split.at(2))), "RMA"));

					} else {

						de_novo_summary_stats_str.join(getAlleleStringAttribute(cur_var->allele(0), "RMA"));
					}

					auto de_novo_summary_stats_emplace = de_novo_summary_stats.emplace(de_novo_summary_stats_str.str(), 0);
					de_novo_summary_stats_emplace.first->second++;
				}
			}

			delete cur_var;

			if ((num_variants % 100000) == 0) {

				std::cout << "[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants" << endl;
			}
		}

		JoiningString variant_summary_stats_header('\t');
		variant_summary_stats_header.join(variant_misc_attributes);

		writeSummaryStats(variant_summary_stats, output_prefix + "_variant.txt", variant_summary_stats_header.str());

		JoiningString mac_summary_stats_header('\t');
		mac_summary_stats_header.join(mac_misc_attributes);

		writeSummaryStats(mac_summary_stats, output_prefix + "_mac.txt", mac_summary_stats_header.str());

		JoiningString allele_summary_stats_header('\t');
		allele_summary_stats_header.join(allele_misc_attributes);

		writeSummaryStats(allele_summary_stats, output_prefix + "_allele.txt", allele_summary_stats_header.str());

		JoiningString de_novo_summary_stats_header('\t');
		de_novo_summary_stats_header.join(de_novo_misc_attributes);

		writeSummaryStats(de_novo_summary_stats, output_prefix + "_de_novo.txt", de_novo_summary_stats_header.str());


		cout << "\n[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants and " << num_alleles << " alleles (including reference and missing)\n" << endl;

		cout << "[" << Utils::getLocalTime() << "] Wrote summary statistics for a total of " << num_variants << " variant, " << num_alleles << " alleles and " << num_genotypes << " genotypes" << endl;
		cout << endl;
	}
}
