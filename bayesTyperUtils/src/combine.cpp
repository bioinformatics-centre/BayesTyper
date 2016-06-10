
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

#include "vcf++/JoiningString.hpp"
#include "vcf++/Auxiliaries.hpp"

#include "combine.hpp"

namespace Combine {

	struct CallSet {

		const string type;
		const string name;
		
		uint num_parsed_variants;
		uint num_excluded_variants;
		uint num_shifted_variants;

		CallSet(const string & type_in, const string & name_in) : type(type_in), name(name_in), num_parsed_variants(0), num_excluded_variants(0), num_shifted_variants(0) {}
	};

    uint leftTrimLength(Allele & allele_1, Allele & allele_2) {

		assert(!(allele_1.isMissing()));
		assert(!(allele_2.isMissing()));

		assert(!(allele_1.seq().empty()));
		assert(!(allele_2.seq().empty()));

		uint left_trim_length = 0;

		if ((allele_1.seq().size() > 1) and (allele_2.seq().size() > 1)) {

			if (allele_1.seq().front() == allele_2.seq().front()) {

				Allele allele_1_tmp = allele_1;
				Allele allele_2_tmp = allele_2;

				left_trim_length = Auxiliaries::partialTrimAllelePair(&allele_1_tmp, &allele_2_tmp).first;
			}
		}

		return left_trim_length;
    }

	void addAnnAlleleAttributes(Variant * cur_var) {

		auto ids = cur_var->ids();
		assert(!(ids.empty()));

		sort(ids.begin(), ids.end());

		JoiningString ids_elements(':');
		ids_elements.join(ids);

		for (uint i = 0; i < cur_var->numAlts(); i++) {

			assert(cur_var->alt(i).info().setValue("AAI", ids_elements.str()));
		} 
	}	

	void addAsmVarAlleleAttributes(Variant * cur_var, const string & call_set_name, const vector<uint> & alt_left_trim_length) {

		assert(alt_left_trim_length.size() == cur_var->numAlts());
		vector<vector<string> > asmvar_asqr_values(cur_var->numAlts());

		for (auto &sid: cur_var->sampleIds()) {

			if (cur_var->getSample(sid).isInformative()) {

				auto genotype = cur_var->getSample(sid).genotypeEstimate();
				assert(genotype.size() == 2);

				assert((genotype.front() == genotype.back()) or (genotype.front() == 0) or (genotype.back() == 0));
				ushort max_allele_idx = max(genotype.front(), genotype.back());

				if (max_allele_idx > 0) {

					auto qr_value = cur_var->getSample(sid).info().getValue<string>("QR");
					assert(qr_value.second);
					assert(qr_value.first != ".");

					auto qr_value_split = Utils::splitString(qr_value.first, '=');
					assert(qr_value_split.size() == 2);

		 			assert(qr_value_split.front().find(",") == string::npos);
		 			assert(qr_value_split.front().find("#") == string::npos);
		 			assert(qr_value_split.front().find(";") == string::npos);

					auto qr_region_split = Utils::splitString(qr_value_split.back(), '-');
					assert(qr_region_split.size() == 2);

					uint region_start_pos = stoi(qr_region_split.front());

					assert(cur_var->allele(max_allele_idx).seq().size() > alt_left_trim_length.at(max_allele_idx - 1));
					region_start_pos += alt_left_trim_length.at(max_allele_idx - 1);

		   		 	JoiningString asqr_value_elements('#');
		    		asqr_value_elements.join(call_set_name);
		    		asqr_value_elements.join(sid);
		    		asqr_value_elements.join(qr_value_split.front());
		    		asqr_value_elements.join(to_string(region_start_pos));	   		 	

					asmvar_asqr_values.at(max_allele_idx - 1).push_back(asqr_value_elements.str());
				}
			}
		}

		for (uint i = 0; i < cur_var->numAlts(); i++) {

			if (!(asmvar_asqr_values.at(i).empty())) {
			
				JoiningString asmvar_asqr_value_elements(':');
				asmvar_asqr_value_elements.join(asmvar_asqr_values.at(i));

				assert(cur_var->alt(i).info().setValue("AsmVar_ASQR", asmvar_asqr_value_elements.str()));

			} else {

				assert(cur_var->alt(i).info().setValue("AsmVar_ASQR", '.'));				
			} 
		}
	}

	vector<uint> getSNPAltIndices(Variant * cur_var) {

		vector<uint> snp_alt_indices;
		snp_alt_indices.reserve(cur_var->numAlts());

		for (uint alt_idx = 0; alt_idx < cur_var->numAlts(); alt_idx++) {

			assert(!(cur_var->alt(alt_idx).isMissing()));

			if (cur_var->alt(alt_idx).seq().size() == cur_var->ref().seq().size()) {

				if (Auxiliaries::alleleAttributes(cur_var->alt(alt_idx), cur_var->ref()).type == Auxiliaries::Type::SNP) {

					snp_alt_indices.push_back(alt_idx);			
				} 
			}
		}

		return snp_alt_indices;
	}

	void updateSortedAtribute(string * out_attribute_str, const string & in_attribute_str) {

		auto out_attribute_str_split = Utils::splitString(*out_attribute_str, ':');
		auto in_attribute_str_split = Utils::splitString(in_attribute_str, ':');

		out_attribute_str_split.reserve(out_attribute_str_split.size() + in_attribute_str_split.size());

		for (auto &in_att: in_attribute_str_split) {

			if (find(out_attribute_str_split.begin(), out_attribute_str_split.end(), in_att) == out_attribute_str_split.end()) {

				out_attribute_str_split.push_back(in_att);
			}
		}

		sort(out_attribute_str_split.begin(), out_attribute_str_split.end());

		JoiningString atribute_elements(':');
		atribute_elements.join(out_attribute_str_split);

		*out_attribute_str = atribute_elements.str();
	}

	void updateAlleleAttribute(Allele & in_allele, Allele * out_allele, const string & allele_info) {
		
		auto in_allele_info_value = in_allele.info().getValue(allele_info);
		assert(in_allele_info_value.second);

		if (in_allele_info_value.first.str() != ".") {

			auto out_allele_info_value = out_allele->info().getValue(allele_info);

			if (out_allele_info_value.second) {

				if (out_allele_info_value.first.str() != ".") {

					auto out_allele_info_value_str = out_allele_info_value.first.str();
					updateSortedAtribute(&out_allele_info_value_str, in_allele_info_value.first.str());
					assert(!(out_allele->info().setValue(allele_info, out_allele_info_value_str)));

				} else {

					assert(!(out_allele->info().setValue(allele_info, in_allele_info_value.first.str())));
				}

			} else {

				assert(out_allele->info().setValue(allele_info, in_allele_info_value.first.str()));
			}
		}
	}

	void updateOriginAttribute(Allele * cur_allele, const string & call_set_name) {

		auto aco_value = cur_allele->info().getValue("ACO");
		
		if (aco_value.second) {

			auto aco_value_str = aco_value.first.str();
			updateSortedAtribute(&aco_value_str, call_set_name);

			assert(!(cur_allele->info().setValue("ACO", aco_value_str)));

		} else {

			assert(cur_allele->info().setValue("ACO", call_set_name));
		}
	}

	void addEmptyAlleleAttribute(Variant * cur_var, const string & allele_info) {

		for (uint i = 0; i < cur_var->numAlts(); i++) {
	
			if (!(cur_var->alt(i).info().hasValue(allele_info))) {

				assert(cur_var->alt(i).info().setValue(allele_info, '.'));
			}
		}
	}

	void appendSeqToAlleles(Variant * cur_var, const string & new_seq) {

	    for (uint allele_idx = 0; allele_idx < cur_var->numAlls(); allele_idx++) {

	        assert(!(cur_var->allele(allele_idx).isMissing()));
	        cur_var->allele(allele_idx).seq().append(new_seq);
	    }
	}

	vector<Variant*> splitAndShiftVariant(Variant * cur_var, const vector<uint> & alt_left_trim_length) {

		assert(cur_var->numSamples() == 0);
		assert(alt_left_trim_length.size() == cur_var->numAlts());

		vector<Variant*> shifted_variants;
		shifted_variants.reserve(cur_var->numAlls());

		vector<uint> shifted_allele_indices;
		shifted_allele_indices.reserve(cur_var->numAlts());

		for (uint alt_idx = 0; alt_idx < cur_var->numAlts(); alt_idx++) {

			if (alt_left_trim_length.at(alt_idx) == 0) {

				continue;
			}
			
			assert(cur_var->ref().seq().front() == cur_var->alt(alt_idx).seq().front());
			assert((cur_var->ref().seq().size() > 1) or (cur_var->alt(alt_idx).seq().size() > 1));

			Allele reference = cur_var->ref();
			Allele alt = cur_var->alt(alt_idx);

			auto trimmed_bases = Auxiliaries::partialTrimAllelePair(&reference, &alt);			

			assert(trimmed_bases.first > 0);
			assert(trimmed_bases.first == alt_left_trim_length.at(alt_idx));

			auto new_var = new Variant(cur_var->chrom(), cur_var->pos() + trimmed_bases.first, reference, vector<Allele>(1, alt), cur_var->info());

			assert(new_var->ids().empty());
			new_var->setIds(cur_var->ids());

			shifted_variants.push_back(new_var);
			shifted_allele_indices.push_back(alt_idx);			
		}

		cur_var->removeAlts(shifted_allele_indices);
		shifted_variants.push_back(cur_var);

		return shifted_variants;
	}

	void addVariant(Variant * cur_var, map<uint, pair<Variant *, set<string> > > * contig_vars, CallSet * call_set, const bool remove_asmvar_snp) {

		call_set->num_parsed_variants++;
		
		vector<uint> snp_alt_indices;
		assert(!(cur_var->ref().seq().empty()));

		if ((call_set->type == "AsmVar") and (remove_asmvar_snp)) {

			snp_alt_indices = getSNPAltIndices(cur_var);
		}

		if (snp_alt_indices.size() == cur_var->numAlts()) {

			call_set->num_excluded_variants++;
			delete cur_var;

		} else {

			cur_var->removeAlts(snp_alt_indices);
			assert(cur_var->numAlts() > 0);

			vector<uint> alt_left_trim_length;
			alt_left_trim_length.reserve(cur_var->numAlts());

			for (uint alt_idx = 0; alt_idx < cur_var->numAlts(); alt_idx++) {

				alt_left_trim_length.push_back(leftTrimLength(cur_var->ref(), cur_var->alt(alt_idx)));
			}

			if (call_set->type == "Ann") {

				addAnnAlleleAttributes(cur_var);
			
			} else if (call_set->type == "AsmVar") {

				addAsmVarAlleleAttributes(cur_var, call_set->name, alt_left_trim_length);
			}

			cur_var->clearSamples();
					
			auto downstream_shifted_vars = splitAndShiftVariant(cur_var, alt_left_trim_length);

			for (auto &shiftet_var: downstream_shifted_vars) {

				if (shiftet_var->numAlts() == 0) {

					assert(shiftet_var->numAlls() == 1);
					delete shiftet_var;

					continue;
				}

				call_set->num_shifted_variants++;

				auto contig_var_it = contig_vars->find(shiftet_var->pos());
				auto out_var = shiftet_var;

				if (contig_var_it == contig_vars->end()) {

					if (call_set->type != "Ann") {

						out_var->setIds(vector<string>());
					}

					out_var->setQual(make_pair(0, false));
					out_var->setFilters(vector<string>());

				} else {

					out_var = contig_var_it->second.first;

					if (call_set->type == "Ann") {

						for (auto &id: shiftet_var->ids()) {

							out_var->addId(id);
						}
					}

					uint min_reference_length = min(out_var->ref().seq().size(), shiftet_var->ref().seq().size());
					assert(out_var->ref().seq().substr(0, min_reference_length) == shiftet_var->ref().seq().substr(0, min_reference_length));

					if (out_var->ref().seq().size() < shiftet_var->ref().seq().size()) {

						appendSeqToAlleles(out_var, shiftet_var->ref().seq().substr(min_reference_length));

					} else if (out_var->ref().seq().size() > shiftet_var->ref().seq().size()) {

						appendSeqToAlleles(shiftet_var, out_var->ref().seq().substr(min_reference_length));
					}

					assert(out_var->ref().seq() == shiftet_var->ref().seq());
				}

				if (contig_var_it == contig_vars->end()) {

					assert(out_var == shiftet_var);

					for (uint alt_idx = 0; alt_idx < shiftet_var->numAlts(); alt_idx++) {

						assert(!(shiftet_var->alt(alt_idx).isMissing()));
						assert(shiftet_var->ref() != shiftet_var->alt(alt_idx));

						updateOriginAttribute(&(shiftet_var->alt(alt_idx)), call_set->name);
					}

					assert(contig_vars->emplace(shiftet_var->pos(), make_pair(shiftet_var, set<string>({call_set->type}))).second);

				} else {

					assert(out_var != shiftet_var);
					assert(out_var == contig_var_it->second.first);

					for (uint alt_idx = 0; alt_idx < shiftet_var->numAlts(); alt_idx++) {

						assert(!(shiftet_var->alt(alt_idx).isMissing()));
						assert(shiftet_var->ref() != shiftet_var->alt(alt_idx));

						auto add_alt_return = out_var->addAlt(shiftet_var->alt(alt_idx));
						updateOriginAttribute(&(add_alt_return.first), call_set->name);

						if (!(add_alt_return.second)) {

							if (call_set->type == "Ann") {

								updateAlleleAttribute(shiftet_var->alt(alt_idx), &(add_alt_return.first), "AAI");
							
							} else if (call_set->type == "AsmVar") {

								updateAlleleAttribute(shiftet_var->alt(alt_idx), &(add_alt_return.first), "AsmVar_ASQR");
							}
						} 
					}

					contig_var_it->second.second.insert(call_set->type);

					delete shiftet_var;
				}
			}
		}
	}
	

	void combine(const vector<string> & in_vcf_files, const string & outfile_prefix, const bool remove_asmvar_snp) {

		cout << "[" << Utils::getLocalTime() << "] Running BayesTyperUtils (" << BTU_VERSION << ") combine on " << in_vcf_files.size() << " files ...\n" << endl;

		assert(in_vcf_files.size() > 0);
		auto num_in_vcfs = in_vcf_files.size();

		vector<pair<unique_ptr<GenotypedVcfFileReader>, CallSet> > in_vcfs;
		in_vcfs.reserve(num_in_vcfs);

		auto in_vcf_file_split = Utils::splitString(in_vcf_files.front(), ':');	
		assert(in_vcf_file_split.size() == 2);

		in_vcfs.emplace_back(unique_ptr<GenotypedVcfFileReader>(new GenotypedVcfFileReader(in_vcf_file_split.back(), true)), CallSet(Utils::splitString(in_vcf_file_split.front(), '_').front(), in_vcf_file_split.front()));

		in_vcfs.back().first->metaData().infoDescriptors().clear();
		Auxiliaries::removeNonRelevantFormatDescriptors(&(in_vcfs.back().first->metaData()), {"GT", "QR"});

		auto output_meta_data = in_vcfs.back().first->metaData();
		auto contigs = output_meta_data.contigs();

		for (uint i = 1; i < num_in_vcfs; i++) {

			in_vcf_file_split = Utils::splitString(in_vcf_files.at(i), ':');
			assert(in_vcf_file_split.size() == 2);

			in_vcfs.emplace_back(unique_ptr<GenotypedVcfFileReader>(new GenotypedVcfFileReader(in_vcf_file_split.back(), true)), CallSet(Utils::splitString(in_vcf_file_split.front(), '_').front(), in_vcf_file_split.front()));

			assert(contigs == in_vcfs.back().first->metaData().contigs());
			
			in_vcfs.back().first->metaData().infoDescriptors().clear();
			Auxiliaries::removeNonRelevantFormatDescriptors(&(in_vcfs.back().first->metaData()), {"GT", "QR"});
		}

		assert(in_vcfs.size() == num_in_vcfs);

		output_meta_data.setFormat("VCFv4.2");
		output_meta_data.clearSamples();
		output_meta_data.miscMeta().clear();
	    output_meta_data.infoDescriptors().clear();
	    output_meta_data.filterDescriptors().clear();
	    output_meta_data.formatDescriptors().clear();

	    vector<pair<string, string> > aco_descriptor_elems({make_pair("ID","ACO"), make_pair("Number","."), make_pair("Type","String"), make_pair("Description","Allele call-set origin (<call-set>:...)")});
	    output_meta_data.infoDescriptors().emplace("ACO", Attribute::DetailedDescriptor(aco_descriptor_elems));

	    unordered_set<string> check_name_uniqueness;

		for (auto &in_vcf: in_vcfs) {

			assert(!(in_vcf.second.type.empty()));
			assert(!(in_vcf.second.name.empty()));

			assert(in_vcf.second.name != ".");
 			assert(in_vcf.second.name.find(",") == string::npos);
 			assert(in_vcf.second.name.find("#") == string::npos);
 			assert(in_vcf.second.name.find(";") == string::npos);			

			assert(check_name_uniqueness.insert(in_vcf.second.name).second);

			if (in_vcf.second.type == "AsmVar") {			

	    		vector<pair<string, string> > asmvar_asqr_descriptor_elems({make_pair("ID","AsmVar_ASQR"), make_pair("Number","."), make_pair("Type","String"), make_pair("Description","AsmVar allele sample query region (<call-set>#<sample>#<contig>#<start_position>:...)")});
			    output_meta_data.infoDescriptors().emplace("AsmVar_ASQR", Attribute::DetailedDescriptor(asmvar_asqr_descriptor_elems));
		 	
		 		for (auto & sample_id: in_vcf.first->metaData().sampleIds()) {

		 			assert(sample_id.find(",") == string::npos);
		 			assert(sample_id.find("#") == string::npos);
		 			assert(sample_id.find(";") == string::npos);
		 		}

			} else if (in_vcf.second.type == "Ann") {

				vector<pair<string, string> > aai_descriptor_elems({make_pair("ID","AAI"), make_pair("Number","A"), make_pair("Type","String"), make_pair("Description","Allele annotation identifiers (<ID>:...)")});
				output_meta_data.infoDescriptors().emplace("AAI", Attribute::DetailedDescriptor(aai_descriptor_elems));	
			}
		}

		VcfFileWriter output_vcf(outfile_prefix + ".vcf", output_meta_data, true);

		unordered_map<string, map<uint, pair<Variant *, set<string> > > > contigs_vars;

		for (auto &contig: contigs) {

			assert(contigs_vars.emplace(contig.id(), map<uint, pair<Variant *, set<string> > >()).second);
		}

		unordered_set<string> combined_contigs_vars;

		uint num_combined_variants = 0;
		vector<pair<uint, uint> > in_vcf_variants_counts(num_in_vcfs, make_pair(0,0));

		for (auto &contig: contigs) {

			for (auto &in_vcf: in_vcfs) {

				Variant * cur_in_var;

				while (in_vcf.first->getNextVariant(&cur_in_var)) {

					if (cur_in_var->chrom() != contig.id()) {

						assert(contigs_vars.count(cur_in_var->chrom()) > 0);
						assert(combined_contigs_vars.count(cur_in_var->chrom()) < 1);

						addVariant(cur_in_var, &(contigs_vars.at(cur_in_var->chrom())), &in_vcf.second, remove_asmvar_snp);
						
						break;
					
					} else {

						addVariant(cur_in_var, &(contigs_vars.at(cur_in_var->chrom())), &in_vcf.second, remove_asmvar_snp);
					}
				}
			}

			assert(combined_contigs_vars.insert(contig.id()).second);

			auto contig_vars_it = contigs_vars.find(contig.id());
			assert(contig_vars_it != contigs_vars.end());

			uint prev_position = 0; 

			auto cit = contig_vars_it->second.begin();

			while (cit != contig_vars_it->second.end()) {

				assert(prev_position < cit->first);
				prev_position = cit->first;

				cit->second.first->removeRedundantAlts();
				Auxiliaries::rightTrimVariant(cit->second.first);

				if (cit->second.second.count("Ann") > 0) {

					addEmptyAlleleAttribute(cit->second.first, "AAI");
				}

				if (cit->second.second.count("AsmVar") > 0) {

					addEmptyAlleleAttribute(cit->second.first, "AsmVar_ASQR");
				}

				num_combined_variants++;
				output_vcf.write(cit->second.first);

				delete cit->second.first;
				cit = contig_vars_it->second.erase(cit);
			}

			cout << "[" << Utils::getLocalTime() << "] Finished chromosome " << contig.id() << endl;
		}

		for (auto &in_vcf: in_vcfs) {

			Variant * cur_in_var;
			assert(!(in_vcf.first->getNextVariant(&cur_in_var)));
		}

		cout << "\n[" << Utils::getLocalTime() << "] Number of parsed variants:\n" << endl;

		for (auto &in_vcf: in_vcfs) {

			cout << "\t- " << in_vcf.second.name << ": " << in_vcf.second.num_parsed_variants << " (" << in_vcf.second.num_excluded_variants << " excluded variants)" << endl;
		}

		cout << "\n[" << Utils::getLocalTime() << "] Number of variants after allele downstream shift (left-trim):\n" << endl;

		for (auto &in_vcf: in_vcfs) {

			cout << "\t- " << in_vcf.second.name << ": " << in_vcf.second.num_shifted_variants << " variants" << endl;
		}

		cout << "\n[" << Utils::getLocalTime() << "] Final number of combined variants: " << num_combined_variants << endl;
		cout << endl;
	}
}
