
/*
annotate.cpp - This file is part of BayesTyper (v0.9)


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


#include <vector>
#include <string>
#include <iostream>
#include <map>
#include <unordered_set>
#include <algorithm>

#include "vcf++/JoiningString.hpp"
#include "vcf++/Auxiliaries.hpp"

#include "annotate.hpp"

namespace Annotate {

	struct AnnotationAllele {

		vector<string> ids;
		Allele alt_allele;
	};

	pair<bool, double> calcReferenceAllelePositionalOverlap(const pair<uint, uint> & ref_allele_1_region, const pair<uint, uint> & ref_allele_2_region) {

		assert(ref_allele_1_region.second >= ref_allele_1_region.first);
		assert(ref_allele_2_region.second >= ref_allele_2_region.first);

		if ((ref_allele_1_region.first == ref_allele_1_region.second) or (ref_allele_2_region.first == ref_allele_2_region.second)) {

			if ((ref_allele_1_region.first == ref_allele_2_region.first) and (ref_allele_1_region.second == ref_allele_2_region.second)) {

				return make_pair(true, 1);
			
			} else {

				return make_pair(false, 0);
			}
		}

		if ((ref_allele_1_region.first == ref_allele_2_region.first) and (ref_allele_1_region.second == ref_allele_2_region.second)) {

			return make_pair(true, 1);
		}

		uint num_base_overlap = 0;

		if (!((ref_allele_1_region.second <= ref_allele_2_region.first) or (ref_allele_2_region.second <= ref_allele_1_region.first))) {

			if (ref_allele_1_region.first <= ref_allele_2_region.first) {

				if (ref_allele_1_region.second <= ref_allele_2_region.second) {

					num_base_overlap = ref_allele_1_region.second - ref_allele_2_region.first;
				
				} else {

					num_base_overlap = ref_allele_2_region.second - ref_allele_2_region.first;				
				}

			} else {

				if (ref_allele_2_region.second <= ref_allele_1_region.second) {

					num_base_overlap = ref_allele_2_region.second - ref_allele_1_region.first;
				
				} else {

					num_base_overlap = ref_allele_1_region.second - ref_allele_1_region.first;					
				}
			}

			assert(num_base_overlap > 0);
		}

		num_base_overlap *= 2;
		
		uint sum_num_nt = ref_allele_1_region.second - ref_allele_1_region.first + ref_allele_2_region.second - ref_allele_2_region.first;
		
		assert(sum_num_nt > 0);
		assert(sum_num_nt > num_base_overlap);

		return make_pair(false, num_base_overlap/static_cast<double>(sum_num_nt));
	}

	uint findMaxSequenceOverlap(const string & shortest_sequence, const string & longest_sequence) {

		assert(shortest_sequence.size() <= longest_sequence.size());

		if (longest_sequence.find(shortest_sequence) != string::npos) {

			return shortest_sequence.size();
		}

		uint max_overlap = 0;

		auto fsit = shortest_sequence.begin();
		fsit++;

		while (fsit != shortest_sequence.end()) {

			if (equal(fsit, shortest_sequence.end(), longest_sequence.begin())) {

				max_overlap = max(max_overlap, static_cast<uint>(shortest_sequence.end() - fsit));
				break;
			}

			fsit++;
		} 

		auto rsit = shortest_sequence.rbegin();
		rsit++;

		while (rsit != shortest_sequence.rend()) {

			if ((shortest_sequence.rend() - rsit) == max_overlap) {

				break;
			}

			if (equal(rsit, shortest_sequence.rend(), longest_sequence.rbegin())) {

				max_overlap = max(max_overlap, static_cast<uint>(shortest_sequence.rend() - rsit));
				break;
			}

			rsit++;
		} 

		return max_overlap;
	}

	pair<bool, double> calcAlternativeAlleleSequenceOverlap(const string & alt_allele_1_seq, const string & alt_allele_2_seq) {
		
		if (alt_allele_1_seq.size() == alt_allele_2_seq.size()) {

			if (alt_allele_1_seq == alt_allele_2_seq) {

				return make_pair(true, 1);
			}
		}

		uint max_nt_overlap = 0;

		if (alt_allele_1_seq.size() <= alt_allele_2_seq.size()) {

			max_nt_overlap = findMaxSequenceOverlap(alt_allele_1_seq, alt_allele_2_seq);

		} else {

			max_nt_overlap = findMaxSequenceOverlap(alt_allele_2_seq, alt_allele_1_seq);
		}

		max_nt_overlap *= 2;

		uint sum_num_nt = alt_allele_1_seq.size() + alt_allele_2_seq.size();
		
		assert(sum_num_nt > 0);
		assert(sum_num_nt > max_nt_overlap);

		return make_pair(false, max_nt_overlap/static_cast<double>(sum_num_nt));
	}

	void annotate(const string & vcf_filename, const string & annotation_filename, const string & output_prefix, const float min_allele_overlap, const bool require_sequence_match) {

		cout << "[" << Utils::getLocalTime() << "] Running BayesTyperUtils (" << BTU_VERSION << ") annotate ...\n" << endl;

		GenotypedVcfFileReader vcf_reader(vcf_filename, true);
		VcfFileReader annotation_reader(annotation_filename, true);

		if (vcf_reader.metaData().infoDescriptors().count("AAI") > 0) {

			assert(vcf_reader.metaData().infoDescriptors().at("AAI").id() == "AAI");
			assert(vcf_reader.metaData().infoDescriptors().at("AAI").number() == Attribute::Number::Dot);
			assert(vcf_reader.metaData().infoDescriptors().at("AAI").type() == Attribute::Type::String);
			assert(vcf_reader.metaData().infoDescriptors().at("AAI").description() == "Allele annotation identifiers (<ID>:...)");

		} else {

	    	vector<pair<string, string> > aai_descriptor_elems({make_pair("ID","AAI"), make_pair("Number","A"), make_pair("Type","String"), make_pair("Description","Allele annotation identifiers (<ID>:...)")});
			assert(vcf_reader.metaData().infoDescriptors().emplace("AAI", Attribute::DetailedDescriptor(aai_descriptor_elems)).second);				
		}

		VcfFileWriter vcf_writer(output_prefix + ".vcf", vcf_reader.metaData(), true);

		auto contigs = vcf_reader.metaData().contigs();
		assert(contigs == annotation_reader.metaData().contigs());

	    annotation_reader.metaData().infoDescriptors().clear();
	    annotation_reader.metaData().filterDescriptors().clear();
	    annotation_reader.metaData().formatDescriptors().clear();

		unordered_map<string, map<pair<uint, uint>, vector<AnnotationAllele> > > annotation_alleles;

		for (auto &contig: contigs) {

			assert(annotation_alleles.emplace(contig.id(), map<pair<uint, uint>, vector<AnnotationAllele> >()).second);
		}

		Variant * cur_vcf_var;
		Variant * cur_anno_var;

		bool cur_vcf_var_not_end = vcf_reader.getNextVariant(&cur_vcf_var);
		bool cur_anno_var_not_end = annotation_reader.getNextVariant(&cur_anno_var);

		uint num_vcf_variants = 0;
		uint num_vcf_alleles = 0;
		uint num_vcf_alleles_missing = 0;

		uint num_annotation_variants = 0;
		uint num_annotation_alleles = 0;

		uint num_annotated_variants = 0;
		uint num_annotated_alleles = 0;

		for (auto &contig: contigs) {

			while (cur_vcf_var_not_end) {

				if (cur_vcf_var->chrom() != contig.id()) {

					break;
				}

				while (cur_anno_var_not_end) {

					if (cur_anno_var->chrom() != contig.id()) {

						break;
					}					

					bool annotation_alleles_buffered = false;

					if ((cur_anno_var->pos() + cur_anno_var->ref().seq().size() - 1) >= cur_vcf_var->pos()) {

						assert(!(cur_anno_var->ids().empty()));
						assert(count(cur_anno_var->ids().begin(), cur_anno_var->ids().end(), ".") == 0);
						assert(count(cur_anno_var->ids().begin(), cur_anno_var->ids().end(), "") == 0);

						auto chrom_annotation_alleles_it = annotation_alleles.find(cur_anno_var->chrom());
						assert(chrom_annotation_alleles_it != annotation_alleles.end());

						for (uint i = 0; i < cur_anno_var->numAlts(); i++) {

							Allele reference_allele = cur_anno_var->ref();
							
							assert(!(cur_anno_var->alt(i).isMissing()));
							
							AnnotationAllele cur_annotation_allele;							
							cur_annotation_allele.ids = cur_anno_var->ids();
							cur_annotation_allele.alt_allele = cur_anno_var->alt(i);

							auto new_pos = cur_anno_var->pos() + Auxiliaries::fullTrimAllelePair(&reference_allele, &cur_annotation_allele.alt_allele).first;

							auto pos_annotation_alleles_it = chrom_annotation_alleles_it->second.emplace(make_pair(new_pos, new_pos + reference_allele.seq().size()), vector<AnnotationAllele>());
							pos_annotation_alleles_it.first->second.push_back(cur_annotation_allele);
						}

						if (cur_anno_var->pos() > (cur_vcf_var->pos() + cur_vcf_var->ref().seq().size() - 1)) {

							annotation_alleles_buffered = true;
						}
					}

					num_annotation_variants++;
					num_annotation_alleles += cur_anno_var->numAlts();

					delete cur_anno_var;
					cur_anno_var_not_end = annotation_reader.getNextVariant(&cur_anno_var);

					if (annotation_alleles_buffered) {

						break;
					}
				}

				while (cur_vcf_var_not_end) {

					if (cur_vcf_var->chrom() != contig.id()) {

						break;
					}

					auto chrom_annotation_alleles_it = annotation_alleles.find(cur_vcf_var->chrom());
					assert(chrom_annotation_alleles_it != annotation_alleles.end());

					if (cur_anno_var_not_end) {

						if (cur_anno_var->chrom() == cur_vcf_var->chrom()) {

							if (!(chrom_annotation_alleles_it->second.empty())) {

								if ((cur_vcf_var->pos() + cur_vcf_var->ref().seq().size() - 1) > chrom_annotation_alleles_it->second.rbegin()->first.first) {

									break;
								}
							}
						}
					}
					
					for (uint i = 0; i < cur_vcf_var->numAlts(); i++) {

						vector<string> aai_values;

						auto cur_vcf_var_aai_val = cur_vcf_var->alt(i).info().getValue<string>("AAI");

						if (cur_vcf_var_aai_val.second) {

							aai_values = Utils::splitString(cur_vcf_var_aai_val.first, ',');

							if (cur_vcf_var->alt(i).isMissing()) {

								assert(aai_values.size() == 1);
								assert(aai_values.front() == ".");
							} 
						}

						if (cur_vcf_var->alt(i).isMissing()) {

							num_vcf_alleles_missing++;
							cur_vcf_var->alt(i).info().setValue<string>("AAI", ".");
							continue;
						}

						pair<Allele, Allele> allele_pair(cur_vcf_var->ref(), cur_vcf_var->alt(i));
						auto new_pos = cur_vcf_var->pos() + Auxiliaries::fullTrimAllelePair(&allele_pair.first, &allele_pair.second).first;

						auto vit = chrom_annotation_alleles_it->second.begin();

						while (vit != chrom_annotation_alleles_it->second.end()) {

							if (cur_vcf_var->pos() > vit->first.second) {

								vit = chrom_annotation_alleles_it->second.erase(vit);

							} else {

								auto reference_allele_overlap = calcReferenceAllelePositionalOverlap(pair<uint, uint>(new_pos, new_pos + allele_pair.first.seq().size()), vit->first);

								if (reference_allele_overlap.first or (reference_allele_overlap.second >= min_allele_overlap)) {

									for (auto &anno_allele: vit->second) {

										pair<bool, double> alternative_allele_overlap = make_pair(true, 1);

										if (require_sequence_match) {

											alternative_allele_overlap = calcAlternativeAlleleSequenceOverlap(allele_pair.second.seq(), anno_allele.alt_allele.seq());
										}

										if (alternative_allele_overlap.first or (alternative_allele_overlap.second >= min_allele_overlap)) {

											for (auto &id: anno_allele.ids) {

												if (reference_allele_overlap.first and alternative_allele_overlap.first) {

													cur_vcf_var->addId(id);
												}

												if (find(aai_values.begin(), aai_values.end(), id) == aai_values.end()) {

													aai_values.push_back(id);
												}
											}
										}
									}
								}

								vit++;
							}
						}

						if (aai_values.empty()) {

							cur_vcf_var->alt(i).info().setValue<string>("AAI", ".");

						} else {

							assert(!(cur_vcf_var->alt(i).isMissing()));
							num_annotated_alleles++;

							sort(aai_values.begin(), aai_values.end());

			   		 		JoiningString aai_value_elements(':');
			    			aai_value_elements.join(aai_values);

							cur_vcf_var->alt(i).info().setValue<string>("AAI", aai_value_elements.str());
						}
					}

					if (!(cur_vcf_var->ids().empty())) {

						num_annotated_variants++;
					}

					vcf_writer.write(cur_vcf_var);

					num_vcf_variants++;
					num_vcf_alleles += cur_vcf_var->numAlts();

					delete cur_vcf_var;
					cur_vcf_var_not_end = vcf_reader.getNextVariant(&cur_vcf_var);
				}	
			}

			while (cur_anno_var_not_end) {

				if (cur_anno_var->chrom() != contig.id()) {

					break;
				}

				num_annotation_variants++;
				num_annotation_alleles += cur_anno_var->numAlts();

				delete cur_anno_var;
				cur_anno_var_not_end = annotation_reader.getNextVariant(&cur_anno_var);
			}

			annotation_alleles.at(contig.id()).clear();

			cout << "[" << Utils::getLocalTime() << "] Finished chromosome " << contig.id() << endl;
		}

		assert(!cur_vcf_var_not_end);
		assert(!(vcf_reader.getNextVariant(&cur_vcf_var)));

		assert(!cur_anno_var_not_end);
		assert(!(annotation_reader.getNextVariant(&cur_anno_var)));

		cout << "\n[" << Utils::getLocalTime() << "] Parsed " << num_vcf_variants << " BayesTyper variants and " << num_vcf_alleles << " alternative alleles (" << num_vcf_alleles_missing << " missing)" << endl;
		cout << "[" << Utils::getLocalTime() << "] Parsed " << num_annotation_variants << " annotation variants and " << num_annotation_alleles << " alternative alleles" << endl;
		cout << "\n[" << Utils::getLocalTime() << "] Annotated " << num_annotated_variants << " BayesTyper variants (perfect match) and " << num_annotated_alleles << " alternative non-missing alleles" << endl;

		cout << endl;
	}
}
