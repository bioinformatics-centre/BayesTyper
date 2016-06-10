
/*
calcFilterStats.cpp - This file is part of BayesTyper (v0.9)


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
#include <algorithm>
#include <regex>

#include "vcf++/JoiningString.hpp"
#include "vcf++/AttributeFilter.hpp"
#include "vcf++/CompareOperators.hpp"
#include "vcf++/ReductionOperator.hpp"
#include "vcf++/SampleAlleleAttributeFilter.hpp"
#include "vcf++/Auxiliaries.hpp"

#include "calcFilterStats.hpp"

namespace CalcFilterStats {

	static const vector<string> stats_header({"VariantType","NumGenotypes", "NumAnnoGenotypes", "NumCalledGenotypes", "NumCalledAnnoGenotypes", "NumVariants", "NumCalledRefAlleles", "NumCalledAltAlleles", "NumCalledAnnoAltAlleles", "InbreedingCount", "InbreedingMean", "InbreedingVar"});

	struct SampleFilters {

		vector<string> ids;
		vector<float> thresholds;

		Attribute::CmpOp * is_less_or_eq;
		Attribute::CmpOp * is_greater_or_eq;
		Attribute::ReductionOp * logical_and_reduction;

		vector<AttributeFilter *> att_filters;
		vector<SampleAlleleAttributeFilter *> sample_att_filters;

		SampleFilters() {

			is_less_or_eq = new Attribute::IsLessOrEqCmpOp();
			is_greater_or_eq = new Attribute::IsGreaterOrEqCmpOp();
			logical_and_reduction = new Attribute::LogicalAndReductionOp();
		}

		void addLessAndEqualFilter(const float threshold, const string & id) {

			ids.push_back(id);
			thresholds.push_back(threshold);

			att_filters.emplace_back(new ValueAttributeFilter(id, Attribute::Value(threshold), is_less_or_eq));
			sample_att_filters.emplace_back(new SampleAlleleAttributeFilter(att_filters.back(), logical_and_reduction));
		}

		void addGreaterAndEqualFilter(const float threshold, const string & id) {

			ids.push_back(id);
			thresholds.push_back(threshold);

			att_filters.emplace_back(new ValueAttributeFilter(id, Attribute::Value(threshold), is_greater_or_eq));
			sample_att_filters.emplace_back(new SampleAlleleAttributeFilter(att_filters.back(), logical_and_reduction));
		}
	};

	struct GenotypeStats {

		uint num_genotypes;
		uint num_anno_genotypes;
		uint num_called_genotypes;
		uint num_called_anno_genotypes;
		uint num_variants;
		uint num_called_ref_allele;
		uint num_called_alt_allele;
		uint num_called_anno_alt_allele;

		uint inbreeding_count;
		double inbreeding_mean;
		double inbreeding_M2;

		GenotypeStats() : num_genotypes(0), num_anno_genotypes(0), num_called_genotypes(0), num_called_anno_genotypes(0), num_variants(0), num_called_ref_allele(0), num_called_alt_allele(0), num_called_anno_alt_allele(0), inbreeding_count(0), inbreeding_mean(0), inbreeding_M2(0) {}
	};


	struct FilterStats {

		SampleFilters sample_filters;
		unordered_map<string, GenotypeStats> variant_type_genotype_stats;

		FilterStats() {

			variant_type_genotype_stats = unordered_map<string, GenotypeStats>({{"SNP", GenotypeStats()}, {"Insertion", GenotypeStats()}, {"Deletion", GenotypeStats()}, {"Inversion", GenotypeStats()}, {"Complex", GenotypeStats()}});
		}
	};


	void calcFilterStats(const string & vcf_filename, const string & output_prefix, const string & parents_trio_regular_expression, const float min_called_probability, const vector<string> & nk_thresholds, const vector<string> & nok_thresholds, const vector<string> & nuk_thresholds) {

		// GenotypedVcfFileReader vcf_reader(vcf_filename, true);
		// Auxiliaries::removeNonRelevantInfoDescriptors(&(vcf_reader.metaData()), {"AAI"});

		// cout << "[" << Utils::getLocalTime() << "] Running BayesTyperUtils (" << BTU_VERSION << ") calcFilterStats on " << vcf_reader.metaData().sampleIds().size() << " samples using " << nk_thresholds.size() * nok_thresholds.size() * nuk_thresholds.size() << " different filter combinations ...\n" << endl;

		// vector<FilterStats> filter_combination_stats;

		// for (auto &thres_nk: nk_thresholds) {

		// 	for (auto &thres_nok: nok_thresholds) {

		// 		for (auto &thres_nuk: nuk_thresholds) {

		// 			filter_combination_stats.emplace_back(FilterStats());
		// 			filter_combination_stats.back().sample_filters.addGreaterAndEqualFilter(stof(thres_nk), "NK");
		// 			filter_combination_stats.back().sample_filters.addGreaterAndEqualFilter(stof(thres_nok), "NOK");
		// 			filter_combination_stats.back().sample_filters.addGreaterAndEqualFilter(stof(thres_nuk), "NUK");
		// 		}
		// 	}
		// }

		// assert(!(filter_combination_stats.empty()));

		// regex parent_sample_id_regex(parents_trio_regular_expression);

		// uint num_parent_samples = 0;

  //   	for (auto &sample_id: vcf_reader.metaData().sampleIds()) {

  //       	if (regex_match(sample_id, parent_sample_id_regex)) {

  //       		num_parent_samples++;
  //       	}
  //       }

		// cout << "[" << Utils::getLocalTime() << "] Calculating population statistics on " << num_parent_samples << " parent samples\n" << endl;

		// Variant * cur_var;

		// uint num_variants = 0;
		// uint num_variants_stats = 0;

		// while (vcf_reader.getNextVariant(&cur_var)) {

		// 	num_variants++;

		// 	bool has_exclusion_attribute = false;

		// 	if (cur_var->numAlls() > 2) {

		// 		has_exclusion_attribute = true;
		// 	}

		// 	assert(cur_var->filters().size());

		// 	if (cur_var->filters().front() != "PASS") {

		// 		has_exclusion_attribute = true;
		// 	}

		// 	if (!has_exclusion_attribute) {

		// 		if (vcf_reader.metaData().getContig(cur_var->chrom()).type() != Contig::Type::Autosomal) {

		// 			has_exclusion_attribute = true;
		// 		}
		// 	}

		// 	if (!has_exclusion_attribute) {

		// 		if (Auxiliaries::hasMissing(*cur_var)) {

		// 			has_exclusion_attribute = true;
		// 		}
		// 	}

		// 	if (!has_exclusion_attribute) {

		// 		assert(cur_var->numAlts() == 1);

		// 		num_variants_stats++;

  //       		for (auto & id: vcf_reader.metaData().sampleIds()) {

  //       			assert((cur_var->getSample(id).ploidy() == Sample::Ploidy::Diploid) and (cur_var->getSample(id).callStatus() == Sample::CallStatus::Complete));
  //       			assert(cur_var->getSample(id).isInformative());
  //       		}

		// 		auto variant_type = Auxiliaries::variantType(*cur_var);
		// 		auto has_annotated_alternative_allele = Auxiliaries::hasAnnotatedAlternativeAllele(*cur_var);

		// 		for (auto &filt_com_stats: filter_combination_stats) {

		// 			GenotypeStats * geno_stats = &(filt_com_stats.variant_type_genotype_stats.at(variant_type));

  //       			for (auto & id: vcf_reader.metaData().sampleIds()) {

  //       				Sample & cur_sample = cur_var->getSample(id);
		// 				bool incl_sample = true;

		//                 for (auto &sample_allele_filter: filt_com_stats.sample_filters.sample_att_filters) {

		//                     if (!sample_allele_filter->pass(cur_sample)) {

		//                         incl_sample = false;
		//                         break;
		//                     }
		//                 }

		//                 if (incl_sample) {

		//                 	auto cur_genotype = cur_sample.genotypeEstimate();
		//                 	assert(cur_genotype.size() == 2);
		//                 	assert((*max_element(cur_genotype.begin(), cur_genotype.end())) <= 1);

		//                 	bool cur_genotype_annotated = true;

		//                 	if (((*max_element(cur_genotype.begin(), cur_genotype.end())) == 1) and !has_annotated_alternative_allele) {

		//                 		cur_genotype_annotated = false;
		//                 	}		         

		//                 	geno_stats->num_genotypes++;

		//                 	if (cur_genotype_annotated) {

		//                 		geno_stats->num_anno_genotypes++;
		//                 	}
	
  //   						auto genotype_gpp = Auxiliaries::getGenotypePosterior(cur_sample);
  //           				assert(genotype_gpp == Auxiliaries::getMaxGenotypePosterior(cur_sample));

  //   						assert(genotype_gpp.second);

		//                 	if (genotype_gpp.first >= min_called_probability) {

		//                 		geno_stats->num_called_genotypes++;

		//                 		if (cur_genotype_annotated) {

		//                 			geno_stats->num_called_anno_genotypes++;
		//                 		}
		//                 	}
		//                 }
		//             }

		// 			auto allele_call_prob = Stats::calcAlleleCallProbAndQualFromAllelePosteriors(cur_var, regex(".+"), filt_com_stats.sample_filters.sample_att_filters).first;

		// 			assert(allele_call_prob.size() == 2);
		// 			assert(cur_var->numAlls() == 2);

		// 			geno_stats->num_variants++;

		// 			if (allele_call_prob.front() >= min_called_probability) {

		// 				geno_stats->num_called_ref_allele++;
		// 			}

		// 			if (allele_call_prob.back() >= min_called_probability) {

		// 				geno_stats->num_called_alt_allele++;

		// 				if (Auxiliaries::isAlleleAnnotated(cur_var->allele(1))) {

		// 					assert(has_annotated_alternative_allele);
		// 					geno_stats->num_called_anno_alt_allele++;
		// 				}
		// 			}

		// 			auto inbreeding_coef = Stats::calcInbreedingCoef(cur_var, parent_sample_id_regex, filt_com_stats.sample_filters.sample_att_filters);

		// 			if (inbreeding_coef.second) {

		// 		        geno_stats->inbreeding_count++;
		// 		        double delta = inbreeding_coef.first - geno_stats->inbreeding_mean;
		// 		        geno_stats->inbreeding_mean += delta/geno_stats->inbreeding_count;
		// 		        geno_stats->inbreeding_M2 += delta * (inbreeding_coef.first - geno_stats->inbreeding_mean);
		// 			}
		// 		}
		// 	}

		// 	delete cur_var;

		// 	if ((num_variants % 100000) == 0) {

		// 		std::cout << "[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants" << endl;
		// 	}
		// }

		// const vector<string> filter_ids = filter_combination_stats.front().sample_filters.ids;

	 //    ofstream stats_writer(output_prefix + "_stats.txt");
	 //    assert(stats_writer.is_open());

		// JoiningString stats_header_elements('\t');
		// stats_header_elements.join(filter_ids);
		// stats_header_elements.join(stats_header);
		// stats_writer << stats_header_elements.str() << endl;

		// for (auto &filt_com_stats: filter_combination_stats) {

		// 	assert(filter_ids == filt_com_stats.sample_filters.ids);
		// 	assert(filter_ids.size() == filt_com_stats.sample_filters.thresholds.size());
		// 	assert(filter_ids.size() == filt_com_stats.sample_filters.att_filters.size());
		// 	assert(filter_ids.size() == filt_com_stats.sample_filters.sample_att_filters.size());

		// 	for (auto &geno_stats: filt_com_stats.variant_type_genotype_stats) {

		// 		JoiningString stats_line_elements('\t');

		// 		for (auto &thres: filt_com_stats.sample_filters.thresholds) {

		// 			stats_line_elements.join(to_string(thres));
		// 		}

		// 		stats_line_elements.join(geno_stats.first);
		// 		stats_line_elements.join(to_string(geno_stats.second.num_genotypes));
		// 		stats_line_elements.join(to_string(geno_stats.second.num_anno_genotypes));
		// 		stats_line_elements.join(to_string(geno_stats.second.num_called_genotypes));
		// 		stats_line_elements.join(to_string(geno_stats.second.num_called_anno_genotypes));
		// 		stats_line_elements.join(to_string(geno_stats.second.num_variants));
		// 		stats_line_elements.join(to_string(geno_stats.second.num_called_ref_allele));
		// 		stats_line_elements.join(to_string(geno_stats.second.num_called_alt_allele));
		// 		stats_line_elements.join(to_string(geno_stats.second.num_called_anno_alt_allele));
		// 		stats_line_elements.join(to_string(geno_stats.second.inbreeding_count));
		// 		stats_line_elements.join(to_string(geno_stats.second.inbreeding_mean));

		// 		if (geno_stats.second.inbreeding_count < 2) {

		// 			stats_line_elements.join("0");

		// 		} else {

		// 			stats_line_elements.join(to_string(geno_stats.second.inbreeding_M2 / (geno_stats.second.inbreeding_count - 1)));
		// 		}

		// 		stats_writer << stats_line_elements.str() << endl;
		// 	}

		// 	delete filt_com_stats.sample_filters.is_less_or_eq;
		// 	delete filt_com_stats.sample_filters.is_greater_or_eq;
		// 	delete filt_com_stats.sample_filters.logical_and_reduction;

		// 	for (auto &filt: filt_com_stats.sample_filters.att_filters) {

		// 		delete filt;
		// 	}

		// 	for (auto &filt: filt_com_stats.sample_filters.sample_att_filters) {

		// 		delete filt;
		// 	}
		// }

		// stats_writer.close();

		// std::cout << "\n[" << Utils::getLocalTime() << "] Completed parsing " << num_variants << " variants" << endl;
		// std::cout << "[" << Utils::getLocalTime() << "] Wrote statistics for " << num_variants_stats << " (remaining excluded)\n" << endl;
	}
}
