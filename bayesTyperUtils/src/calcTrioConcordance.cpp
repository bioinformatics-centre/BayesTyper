
/*
calcTrioConcordance.cpp - This file is part of BayesTyper (v0.9)


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

#include "calcTrioConcordance.hpp"

namespace CalcTrioConcordance {

	struct ConcordantCounts {

		ulong count;
		ulong concordant;

		ConcordantCounts() {

			count = 0;
			concordant = 0;
		}
	};

	void writeTypeConcordance(const map<string, ConcordantCounts> & type_concordance_counts) {

		for (auto & type : type_concordance_counts) {

			cout << "\t\t- " << type.first << ": " << static_cast<float>(type.second.concordant)/(type.second.count) << " (" << type.second.count << ")\n";
		}
	}

	void calcTrioConcordance(const string & vcf_filename, const string & output_prefix, const string & trio_info_str) {

        cout << "[" << Utils::getLocalTime() << "] Running BayesTyperUtils (" << BTU_VERSION << ") calcTrioConcordance ...\n" << endl;

		GenotypedVcfFileReader vcf_reader(vcf_filename, true);	

		Auxiliaries::removeNonRelevantInfoDescriptors(&(vcf_reader.metaData()), {"ACO", "RMA"});
		Auxiliaries::removeNonRelevantFormatDescriptors(&(vcf_reader.metaData()), {"GT", "GPP"});

		vector<Trio::TrioInfo> all_trio_info;

		if (trio_info_str.empty()) {

			all_trio_info = Trio::parseGenomeDKPedigree(vcf_reader.metaData());

		} else {

			all_trio_info = Trio::parsePedigree(vcf_reader.metaData(), trio_info_str);
		}

		map<string, ConcordantCounts> trio_concordance;

		cout << "[" << Utils::getLocalTime() << "] Assessing concordance between:\n" << endl;

		for (auto &cur_trio_info: all_trio_info) {

			cout << "\t- " << cur_trio_info.id << ": Father " << cur_trio_info.father << ", mother " << cur_trio_info.mother << ", and child " << cur_trio_info.child << endl;

			assert(trio_concordance.emplace(cur_trio_info.id, ConcordantCounts()).second);
		}

		cout << endl;

		ulong num_variants = 0;
		ulong num_skipped_chr = 0;
		ulong num_filtered_trios = 0;
		ulong num_uninformative_trios = 0;

		map<string, ConcordantCounts> variant_type_concordance;
		map<string, ConcordantCounts> chromosome_type_concordance({{"Autosomal", ConcordantCounts()}, {"ChrX", ConcordantCounts()}, {"ChrY", ConcordantCounts()}});

		unordered_map<string, pair<ulong, ulong> > concordance_trio_stats;

		Variant * cur_var;

		while (vcf_reader.getNextVariant(&cur_var)) {

			num_variants++;

			if ((vcf_reader.metaData().getContig(cur_var->chrom()).type() != Contig::Type::Autosomal) and (vcf_reader.metaData().getContig(cur_var->chrom()).type() != Contig::Type::ChrX) and (vcf_reader.metaData().getContig(cur_var->chrom()).type() != Contig::Type::ChrY)) {

				num_skipped_chr++;
			
			} else {

				const string chromosome_type = vcf_reader.metaData().getContig(cur_var->chrom()).typeStr();
				const string variant_type = Auxiliaries::variantType(*cur_var);
				
				const bool has_missing = Auxiliaries::hasMissing(*cur_var);
				const bool has_repeat = Auxiliaries::hasRepeat(*cur_var);

	        	uint max_allele_length = 0;
	            uint max_num_ambiguous = 0;
	        	uint max_abs_sv_length = 0;

				for (uint alt_idx = 0; alt_idx < cur_var->numAlts(); alt_idx++) {

					auto allele_attributes = Auxiliaries::alleleAttributes(cur_var->alt(alt_idx), cur_var->ref());
				
					max_allele_length = max(max_allele_length, allele_attributes.length);
					max_num_ambiguous = max(max_num_ambiguous, allele_attributes.num_ambiguous);
					max_abs_sv_length = max(max_abs_sv_length, static_cast<uint>(abs(allele_attributes.sv_length)));
				}

				assert(cur_var->filters().size() == 1);
				
				if (variant_type_concordance.count(variant_type) < 1) {

					assert(variant_type_concordance.emplace(variant_type, ConcordantCounts()).second);
				}

				string variant_origin = "Unknown";

				assert(cur_var->numAlts() > Auxiliaries::hasMissing(*cur_var));
        		
    			if (cur_var->numAlts() == (1 + static_cast<uint>(Auxiliaries::hasMissing(*cur_var)))) {

					auto aco_value = cur_var->alt(0).info().getValue<string>("ACO");
					
					if (aco_value.second) {

						variant_origin = aco_value.first; 			
					}
    			}

				for (auto &cur_trio_info: all_trio_info) {

					Trio cur_trio = Trio(*cur_var, cur_trio_info);

					if (cur_trio.isFiltered()) {

						num_filtered_trios++;

					} else if (!cur_trio.isInformative()) {

						num_uninformative_trios++;

					} else {

						trio_concordance.at(cur_trio_info.id).count++;
						variant_type_concordance.at(variant_type).count++;
						chromosome_type_concordance.at(chromosome_type).count++;

						if (cur_trio.isConcordant()) {

							trio_concordance.at(cur_trio_info.id).concordant++;
							variant_type_concordance.at(variant_type).concordant++;
							chromosome_type_concordance.at(chromosome_type).concordant++;
						}

						JoiningString trio_stats_line_elements('\t');
						trio_stats_line_elements.join(chromosome_type);
						trio_stats_line_elements.join(cur_var->filters().front());
						trio_stats_line_elements.join(cur_trio_info.id);
						trio_stats_line_elements.join(Utils::floatToString(cur_trio.minGPP(), 2));
						trio_stats_line_elements.join(variant_type);
						trio_stats_line_elements.join(variant_origin);
						trio_stats_line_elements.join(Utils::boolToString(has_missing));
						trio_stats_line_elements.join(Utils::boolToString(has_repeat));
						trio_stats_line_elements.join(Utils::boolToString(cur_trio.isReferenceCall()));
						trio_stats_line_elements.join(Utils::boolToString(cur_trio.isParentsBiAllelicHeterozygote()));
						trio_stats_line_elements.join(to_string(cur_var->numAlls()));
						trio_stats_line_elements.join(to_string(max_allele_length));
						trio_stats_line_elements.join(to_string(max_num_ambiguous));
						trio_stats_line_elements.join(to_string(max_abs_sv_length));

						auto concordance_trio_stats_emplace = concordance_trio_stats.emplace(trio_stats_line_elements.str(), pair<ulong, ulong>(0,0));

						concordance_trio_stats_emplace.first->second.first++;
						concordance_trio_stats_emplace.first->second.second += cur_trio.isConcordant();
					}
				}
			}
			
			delete cur_var;

			if ((num_variants % 100000) == 0) {

				std::cout << "[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants" << endl;
			}
		}

	    ofstream concordance_trio_stats_writer(output_prefix + "_concordance_trio_stats.txt");
	    assert(concordance_trio_stats_writer.is_open());

		concordance_trio_stats_writer << "NumGenotypes\tNumConcordantGenotypes\tChromosomeType\tFilter\tTrioId\tMinGPP\tVariantType\tOrigin\tHasMissing\tHasRepeat\tIsReferenceCall\tIsParentsBiAllelicHeterozygote\tNumAlleles\tMaxAlleleLength\tMaxNumAmbiguous\tMaxAbsAlleleSVLength\n";

		for (auto &stats: concordance_trio_stats) {

			concordance_trio_stats_writer << stats.second.first << "\t" << stats.second.second << "\t" << stats.first << "\n";
		}

		concordance_trio_stats_writer.close();

		cout << "\n[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants:\n" << endl;
		cout << "\t\t- Skipped " << num_skipped_chr << " variants on non autosomal or allosomal chromosome" << endl;

		ulong remaining_variants = num_variants - num_skipped_chr;

		cout << "\n[" << Utils::getLocalTime() << "] Parsed " << remaining_variants * all_trio_info.size() << " trio genotypes:\n" << endl;
		cout << "\t\t- Skipped " << num_filtered_trios << " filtered trio genotypes (at least one genotype filtered)" << endl;
		cout << "\t\t- Skipped " << num_uninformative_trios << " uninformative trio genotypes (trios with female child on chromosome Y or trios with incorrect ploidy)" << endl;

		ulong num_trios = 0;
		ulong num_concordant = 0;

		for (auto & trio : trio_concordance) {

			num_trios += trio.second.count;
			num_concordant += trio.second.concordant;
		}

		assert((remaining_variants * all_trio_info.size() - num_filtered_trios - num_uninformative_trios) == num_trios);

		cout << "\n\tOverall trio concordance:\n" << endl;
		cout << "\t\t- Sensitivity: " << static_cast<float>(num_trios)/(remaining_variants * all_trio_info.size() - num_uninformative_trios) << endl;
		cout << "\t\t- Accuracy: " << static_cast<float>(num_concordant)/(num_trios) << endl;

		cout << "\n\tTrio concordance:\n" << endl;
		writeTypeConcordance(trio_concordance);

		cout << "\n\tVariant type trio concordance:\n" << endl;
		writeTypeConcordance(variant_type_concordance);

		cout << "\n\tChromosome type trio concordance:\n" << endl;
		writeTypeConcordance(chromosome_type_concordance);

		cout << "\n" << endl;
	}
}
