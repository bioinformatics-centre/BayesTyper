
/*
combine.cpp - This file is part of BayesTyper (v0.9)


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
#include <unordered_map>
#include <set>
#include <sstream>
#include <iostream>
#include <algorithm>

#include "JoiningString.hpp"
#include "Auxiliaries.hpp"

#include "combine.hpp"

namespace Combine {

	struct CallSet {

		const string name;

		VcfFileReader * vcf_reader;
		Variant * variant;

		uint num_alt_alleles;

		CallSet(const string & name_in, const string & vcf_filename) : name(name_in), vcf_reader(nullptr), variant(nullptr), num_alt_alleles(0) {

			vcf_reader = new VcfFileReader(vcf_filename, true);
		}
	};

	struct AlleleSet {

		vector<pair<Variant *, uint> > alleles;
		uint prev_end_pos;
		
		string alt_seq;
		string ref_seq;
	};

	string mergeOriginValues(const string & aco_values_1, const string & aco_values_2) {

		assert(!(aco_values_1.empty()));
		assert(!(aco_values_2.empty()));

		assert(aco_values_1 != ".");
		assert(aco_values_2 != ".");

		vector<string> aco_value_1_split = Utils::splitString(aco_values_1, ':');
		vector<string> aco_value_2_split = Utils::splitString(aco_values_2, ':');

		for (auto & aco_value_2: aco_value_2_split) {

			if (find(aco_value_1_split.begin(), aco_value_1_split.end(), aco_value_2) == aco_value_1_split.end()) {

				aco_value_1_split.push_back(aco_value_2);
			}
		}

		sort(aco_value_1_split.begin(), aco_value_1_split.end());

		JoiningString merged_aco_values(':');
		merged_aco_values.join(aco_value_1_split);

		return merged_aco_values.str();
	}

	void getRedundantAlleleSets(vector<AlleleSet> * redundant_allele_sets, const AlleleSet & allele_set, const string & ref_allele_seq, const string & alt_allele_seq, typename map<uint, Variant *>::iterator contig_variants_it, typename map<uint, Variant *>::iterator contig_variants_eit) {

		if (contig_variants_it != contig_variants_eit) {

			string inter_var_ref_seq_1 = "";
			string inter_var_ref_seq_2 = "";

	    	if (!(allele_set.alleles.empty()) and (allele_set.prev_end_pos < contig_variants_it->second->pos())) {

	    		inter_var_ref_seq_1 = ref_allele_seq.substr(allele_set.ref_seq.size(), contig_variants_it->second->pos() - allele_set.prev_end_pos - 1);
	    		inter_var_ref_seq_2 = alt_allele_seq.substr(allele_set.alt_seq.size(), contig_variants_it->second->pos() - allele_set.prev_end_pos - 1);
	    	}

	    	if (inter_var_ref_seq_1 == inter_var_ref_seq_2) {

				for (uint allele_idx = 0; allele_idx < contig_variants_it->second->numAlls(); allele_idx++) {

			    	auto cur_allele_set = allele_set;
			    	cur_allele_set.alleles.emplace_back(contig_variants_it->second, allele_idx);

			    	Allele cur_ref_allele = contig_variants_it->second->ref();
			    	Allele cur_allele = contig_variants_it->second->allele(allele_idx);

			    	assert(!(cur_ref_allele.isMissing()));
			    	assert(!(cur_allele.isMissing()));

			    	Auxiliaries::rightTrimAllelePair(&cur_ref_allele, &cur_allele);

			    	assert(!(cur_ref_allele.seq().empty()));
			    	assert(!(cur_allele.seq().empty()));	

			    	cur_allele_set.prev_end_pos = contig_variants_it->second->pos() + cur_ref_allele.seq().size() - 1;

	    			cur_allele_set.ref_seq.append(inter_var_ref_seq_1);
	    			cur_allele_set.alt_seq.append(inter_var_ref_seq_1);

			    	cur_allele_set.ref_seq.append(cur_ref_allele.seq());
			    	cur_allele_set.alt_seq.append(cur_allele.seq());

			    	if ((cur_allele_set.ref_seq == ref_allele_seq) and (cur_allele_set.alt_seq == alt_allele_seq)) {

			    		if (cur_allele_set.alleles.size() > 1) {
			    			
			    			redundant_allele_sets->push_back(cur_allele_set);	
			    		}		    			
			    	
			    	} else if ((cur_allele_set.ref_seq == ref_allele_seq.substr(0, cur_allele_set.ref_seq.size())) and (cur_allele_set.alt_seq == alt_allele_seq.substr(0, cur_allele_set.alt_seq.size()))) {

			    		auto next_contig_variants_it = contig_variants_it;
			    		next_contig_variants_it++;

			    		getRedundantAlleleSets(redundant_allele_sets, cur_allele_set, ref_allele_seq, alt_allele_seq, next_contig_variants_it, contig_variants_eit);
			    	}
				}	
			}
		}
	}

    bool isAltAlleleRedundant(typename map<uint, Variant *>::iterator contig_variants_it, typename map<uint, Variant *>::iterator contig_variants_eit, const uint alt_allele_idx) {

    	assert(contig_variants_it != contig_variants_eit);
    	
    	Allele ref_allele = contig_variants_it->second->ref();
    	Allele alt_allele = contig_variants_it->second->alt(alt_allele_idx);

    	assert(!(ref_allele.isMissing()));
    	assert(!(alt_allele.isMissing()));

    	Auxiliaries::rightTrimAllelePair(&ref_allele, &alt_allele);

    	assert(!(ref_allele.seq().empty()));
    	assert(!(alt_allele.seq().empty()));

    	if (ref_allele.seq().size() > 1) {

    		AlleleSet allele_set;
    		vector<AlleleSet> redundant_allele_sets;

    		getRedundantAlleleSets(&(redundant_allele_sets), allele_set, ref_allele.seq(), alt_allele.seq(), contig_variants_it, contig_variants_eit);

			if (!(redundant_allele_sets.empty())) {

				auto aco_value_original_allele = alt_allele.info().getValue<string>("ACO");
				assert(aco_value_original_allele.second);

				for (auto & redundant_allele_set: redundant_allele_sets) {

					assert(redundant_allele_set.alleles.size() > 1);
	    			assert(redundant_allele_set.ref_seq == ref_allele.seq());
	    			assert(redundant_allele_set.alt_seq == alt_allele.seq());

					for (auto & redundant_allele: redundant_allele_set.alleles) {

						if (redundant_allele.second > 0) {

							auto aco_value_redundant_allele = redundant_allele.first->allele(redundant_allele.second).info().getValue<string>("ACO");
							assert(aco_value_redundant_allele.second);

							redundant_allele.first->allele(redundant_allele.second).info().setValue<string>("ACO", mergeOriginValues(aco_value_redundant_allele.first, aco_value_original_allele.first));
						}
					}
				}

				return true;
			}
    	}

    	return false;
	}

	void updateOriginAttribute(Allele * cur_allele, const string & call_set_name) {

		assert(!(cur_allele->isMissing()));
		auto aco_value = cur_allele->info().getValue<string>("ACO");

		if (aco_value.second) {
	
			auto aco_value_split = Utils::splitString(aco_value.first, ':');

			if (find(aco_value_split.begin(), aco_value_split.end(), call_set_name) == aco_value_split.end()) {

				aco_value_split.push_back(call_set_name);
			}

			sort(aco_value_split.begin(), aco_value_split.end());

			JoiningString aco_value_elements(':');
			aco_value_elements.join(aco_value_split);

			assert(!(cur_allele->info().setValue<string>("ACO", aco_value_elements.str())));

		} else {

			assert(cur_allele->info().setValue<string>("ACO", call_set_name));
		}
	}

	void appendSeqToAlleles(Variant * cur_var, const string & new_seq) {

	    for (uint allele_idx = 0; allele_idx < cur_var->numAlls(); allele_idx++) {

	        assert(!(cur_var->allele(allele_idx).isMissing()));
	        cur_var->allele(allele_idx).seq().append(new_seq);
	    }
	}

	uint addVariant(Variant * cur_var, map<uint, Variant *> * contig_variants, const string & call_set_name, const bool exclude_ambiguous_alleles) {

		assert(!(cur_var->ref().seq().empty()));
		assert(cur_var->numAlts() > 0);

		vector<uint> excluded_alt_allele_indices;
		excluded_alt_allele_indices.reserve(cur_var->numAlts());

        if (exclude_ambiguous_alleles and (cur_var->ref().seq().find_first_not_of("ACGT") != string::npos)) {

        	excluded_alt_allele_indices = vector<uint>(cur_var->numAlts());
        	iota(excluded_alt_allele_indices.begin(), excluded_alt_allele_indices.end(), 0);
        } 

        assert(excluded_alt_allele_indices.empty() or (excluded_alt_allele_indices.size() == cur_var->numAlts()));

        if (excluded_alt_allele_indices.empty()) {

			for (uint alt_allele_idx = 0; alt_allele_idx < cur_var->numAlts(); alt_allele_idx++) {

				if (cur_var->alt(alt_allele_idx).isID()) {

					excluded_alt_allele_indices.push_back(alt_allele_idx);

				} else if (cur_var->alt(alt_allele_idx).isMissing()) {

					excluded_alt_allele_indices.push_back(alt_allele_idx);
				
				} else if (exclude_ambiguous_alleles and (cur_var->alt(alt_allele_idx).seq().find_first_not_of("ACGT") != string::npos)) {

					excluded_alt_allele_indices.push_back(alt_allele_idx);
				}
			}
        }

		assert(excluded_alt_allele_indices.size() <= cur_var->numAlts());
		cur_var->removeAlts(excluded_alt_allele_indices);

		const uint num_alt_alleles = cur_var->numAlts();

		auto contig_variants_it = contig_variants->emplace(cur_var->pos(), cur_var);

		if (contig_variants_it.second) {

			contig_variants_it.first->second->setIds(vector<string>());
			contig_variants_it.first->second->setQual(make_pair(0, false));
			contig_variants_it.first->second->setFilters(vector<string>());

			for (uint alt_allele_idx = 0; alt_allele_idx < contig_variants_it.first->second->numAlts(); alt_allele_idx++) {

				updateOriginAttribute(&(contig_variants_it.first->second->alt(alt_allele_idx)), call_set_name);
			}

		} else {

			uint min_ref_length = min(contig_variants_it.first->second->ref().seq().size(), cur_var->ref().seq().size());
			assert(contig_variants_it.first->second->ref().seq().substr(0, min_ref_length) == cur_var->ref().seq().substr(0, min_ref_length));

			if (contig_variants_it.first->second->ref().seq().size() < cur_var->ref().seq().size()) {

				appendSeqToAlleles(contig_variants_it.first->second, cur_var->ref().seq().substr(min_ref_length));

			} else if (contig_variants_it.first->second->ref().seq().size() > cur_var->ref().seq().size()) {

				appendSeqToAlleles(cur_var, contig_variants_it.first->second->ref().seq().substr(min_ref_length));
			}

			assert(contig_variants_it.first->second->ref().seq() == cur_var->ref().seq());

			for (uint alt_allele_idx = 0; alt_allele_idx < cur_var->numAlts(); alt_allele_idx++) {

				auto alt_allele_it = contig_variants_it.first->second->addAlt(cur_var->alt(alt_allele_idx));
				updateOriginAttribute(&(alt_allele_it.first), call_set_name);
			}

			delete cur_var;
		}

		return num_alt_alleles;
	}

	void combine(const vector<string> & in_vcf_files, const string & outfile_prefix, const bool exclude_ambiguous_alleles) {

		cout << "[" << Utils::getLocalTime() << "] Running BayesTyperTools (" << BT_VERSION << ") combine on " << in_vcf_files.size() << " files ...\n" << endl;

		assert(!(in_vcf_files.empty()));

		vector<CallSet> input_call_sets;
		input_call_sets.reserve(in_vcf_files.size());

		for (auto & in_vcf_file: in_vcf_files) {

			auto in_vcf_file_split = Utils::splitString(in_vcf_file, ':');
			assert(in_vcf_file_split.size() == 2);

			input_call_sets.emplace_back(in_vcf_file_split.front(), in_vcf_file_split.back());
		}

		assert(!(input_call_sets.empty()));
		assert(input_call_sets.size() == in_vcf_files.size());

		for (auto & input_call_set: input_call_sets) {

			input_call_set.vcf_reader->getNextVariant(&(input_call_set.variant));

			assert(!(input_call_set.name.empty()));
			assert(input_call_set.name != ".");
 			assert(input_call_set.name.find(",") == string::npos);
 			assert(input_call_set.name.find("#") == string::npos);
 			assert(input_call_set.name.find(";") == string::npos);

			input_call_set.vcf_reader->metaData().infoDescriptors().clear();
			input_call_set.vcf_reader->metaData().formatDescriptors().clear();

			assert(input_call_set.vcf_reader->metaData().contigs() == input_call_sets.front().vcf_reader->metaData().contigs());
		}

		auto output_meta_data = input_call_sets.front().vcf_reader->metaData();

		output_meta_data.setFormat("VCFv4.2");
		output_meta_data.miscMeta().clear();
	    output_meta_data.infoDescriptors().clear();
	    output_meta_data.filterDescriptors().clear();
	    output_meta_data.formatDescriptors().clear();
	    output_meta_data.clearSamples();

	    vector<pair<string, string> > aco_descriptor({make_pair("ID","ACO"), make_pair("Number","A"), make_pair("Type","String"), make_pair("Description","Allele call set origin (<call set>:...)")});
	    output_meta_data.infoDescriptors().emplace("ACO", Attribute::DetailedDescriptor(aco_descriptor));

		VcfFileWriter vcf_writer(outfile_prefix + ".vcf", output_meta_data, true);

		auto contigs = output_meta_data.contigs();
		uint num_combined_alt_alleles = 0;

		for (auto &contig: contigs) {

			map<uint, Variant *> contig_variants;

			for (auto & input_call_set: input_call_sets) {

				while (input_call_set.variant) {
					
					if (input_call_set.variant->chrom() != contig.id()) {

						break;
					}

					input_call_set.num_alt_alleles += addVariant(input_call_set.variant, &contig_variants, input_call_set.name, exclude_ambiguous_alleles);
					input_call_set.vcf_reader->getNextVariant(&(input_call_set.variant));
				}
			}

			auto contig_variants_it = contig_variants.begin();

			while (contig_variants_it != contig_variants.end()) {

				assert(contig_variants_it->first == contig_variants_it->second->pos());

				vector<uint> redundant_alt_allele_indices;

				for (uint alt_allele_idx = 0; alt_allele_idx < contig_variants_it->second->numAlts(); alt_allele_idx++) {

					assert(!(contig_variants_it->second->alt(alt_allele_idx).isID()));
					assert(!(contig_variants_it->second->alt(alt_allele_idx).isMissing()));
					assert(contig_variants_it->second->alt(alt_allele_idx) != contig_variants_it->second->ref());

            		if (isAltAlleleRedundant(contig_variants_it, contig_variants.end(), alt_allele_idx)) {

            			redundant_alt_allele_indices.push_back(alt_allele_idx);                			
            		}
				}

				assert(redundant_alt_allele_indices.size() <= contig_variants_it->second->numAlts());
				contig_variants_it->second->removeAlts(redundant_alt_allele_indices);
				
				if (contig_variants_it->second->numAlts() > 0) {
					
					Auxiliaries::rightTrimVariant(contig_variants_it->second);

					num_combined_alt_alleles += contig_variants_it->second->numAlts();
					vcf_writer.write(contig_variants_it->second);
				}

				delete contig_variants_it->second;
				contig_variants_it++;
			}

			cout << "[" << Utils::getLocalTime() << "] Finished chromosome " << contig.id() << endl;
		}

		cout << "\n[" << Utils::getLocalTime() << "] Number of included alternative alleles (id";

		if (exclude_ambiguous_alleles) {

			cout << ", missing and ambiguous";
		
		} else {

			cout << " and missing";
		}

		cout << "):\n" << endl;

		for (auto & input_call_set: input_call_sets) {

			assert(!(input_call_set.vcf_reader->getNextVariant(&(input_call_set.variant))));
			assert(!(input_call_set.variant));

			delete input_call_set.vcf_reader;

			cout << "\t- " << input_call_set.name << ": " << input_call_set.num_alt_alleles << endl;
		}

		cout << "\n[" << Utils::getLocalTime() << "] Number of alternative alleles in the combined set: " << num_combined_alt_alleles << endl;
		cout << endl;
	}
}
