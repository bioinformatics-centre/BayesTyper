
/*
addAttributes.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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
#include <unordered_map>
#include <regex>

#include "VcfFile.hpp"
#include "Variant.hpp"
#include "JoiningString.hpp"
#include "Trio.hpp"
#include "Auxiliaries.hpp"
#include "FastaRecord.hpp"
#include "FastaReader.hpp"
#include "Stats.hpp"

namespace AddAttributes {

    struct SortCover {

        bool operator() (const pair<uint, uint> & lhs, const pair<uint, uint> & rhs) {

            return lhs.first < rhs.first;
        }
    };

    unordered_map<string, map<string, vector<pair<uint, uint> > > > parseRepeatFile(const string & repeat_filename) {

        unordered_map<string, map<string, vector<pair<uint, uint> > > > repeat_annotations;

        ifstream repeat_file(repeat_filename);

        string line;
        uint line_count = 0;

        while (getline(repeat_file, line)) {

            line_count++;

            if (line_count > 3) {

                auto line_split = Utils::splitStringEmptyIgnore(line, ' ');
                assert((line_split.size() == 15) or (line_split.size() == 16));

                string name = line_split.at(10);
                assert(!(name.empty()));

                assert(name.find(":") == string::npos);
                assert(name.find("#") == string::npos);

                assert(line_split.at(5).front() != '(');
                assert(line_split.at(5).back() != ')');
                assert(line_split.at(6).front() != '(');
                assert(line_split.at(6).back() != ')');
                assert(line_split.at(7).front() == '(');
                assert(line_split.at(7).back() == ')');

                uint query_start = stoi(line_split.at(5));
                uint query_end = stoi(line_split.at(6));

                assert(query_start <= query_end);

                auto repeat_annotations_it = repeat_annotations.emplace(line_split.at(4), map<string, vector<pair<uint, uint> > >());
                
                auto repeat_annotations_family_it = repeat_annotations_it.first->second.emplace(name, vector<pair<uint, uint> >());
                repeat_annotations_family_it.first->second.emplace_back(query_start, query_end);

                auto repeat_annotations_all_it = repeat_annotations_it.first->second.emplace("All", vector<pair<uint, uint> >());
                repeat_annotations_all_it.first->second.emplace_back(query_start, query_end);
            }
        }

        repeat_file.close();

        return repeat_annotations;
    }

    uint calculateNucleotideCover(vector<pair<uint, uint> > cover) {

        assert(!(cover.empty()));
        sort(cover.begin(), cover.end(), SortCover());

        auto cover_it = cover.begin();
        assert(cover_it->first <= cover_it->second);

        uint nt_cover = cover_it->second - cover_it->first + 1;
        uint prev_end = cover_it->second;

        cover_it++;

        while (cover_it != cover.end()) {

            assert(cover_it->first <= cover_it->second);

            if (prev_end < cover_it->second) {

                if (prev_end < cover_it->first) {

                    nt_cover += cover_it->second - cover_it->first + 1;

                } else {

                    nt_cover += cover_it->second - prev_end;
                }

                prev_end = cover_it->second;
            }

            cover_it++;
        }

        return nt_cover;
    }

	void addAttributes(const string & vcf_filename, const string & output_prefix, const string & genome_filename, const string & repeat_filename, const string & indepedent_samples_regex_str, const string & trio_sample_info_str) {


		cout << "[" << Utils::getLocalTime() << "] Running BayesTyperTools (" << BT_VERSION << ") addAttributes ...\n" << endl;

        unordered_map<string, FastaRecord*> genome_seqs;
        unordered_map<string, map<string, vector<pair<uint, uint> > > > repeat_annotations;
        regex indepedent_samples_regex;
		vector<Trio::TrioInfo> all_trio_info;

		GenotypedVcfFileReader vcf_reader(vcf_filename, true);
		auto output_meta_data = vcf_reader.metaData();

		const vector<string> sample_ids = output_meta_data.sampleIds();

		cout << "[" << Utils::getLocalTime() << "] Adding the following attributes:\n" << endl;

		assert(!(genome_filename.empty()) or !(repeat_filename.empty()) or !(indepedent_samples_regex_str.empty()) or !(trio_sample_info_str.empty()));

		if (!(genome_filename.empty())) {

			cout << "\t - Homopolymer length (HPL)" << endl;

	        FastaReader genome_reader(genome_filename);
	        FastaRecord * cur_fasta_rec;

	        while (genome_reader.getNextRecord(&cur_fasta_rec)) {

	            assert(genome_seqs.emplace(cur_fasta_rec->id(), cur_fasta_rec).second);
	            genome_seqs.at(cur_fasta_rec->id())->convertToUppercase();
	        }

			output_meta_data.infoDescriptors().emplace("HPL", Attribute::DetailedDescriptor("HPL", Attribute::Number::One, Attribute::Type::String, "Homopolymer length (<Length>:<Nucleotide>)"));
		}

		if (!(repeat_filename.empty())) {

			cout << "\t - RepeatMasker annotations (RMA)" << endl;

			repeat_annotations = parseRepeatFile(repeat_filename);

			output_meta_data.infoDescriptors().emplace("RMA", Attribute::DetailedDescriptor("RMA", Attribute::Number::A, Attribute::Type::String, "RepeatMasker annotations (<family#nucleotide_cover>:...)"));
		}

		if (!(indepedent_samples_regex_str.empty())) {

			cout << "\t - Absolute inbreeding coefficient (IBC)" << endl;

			indepedent_samples_regex = regex(indepedent_samples_regex_str);

			output_meta_data.infoDescriptors().emplace("IBC", Attribute::DetailedDescriptor("IBC", Attribute::Number::One, Attribute::Type::String, "Absolute inbreeding coefficient (<Coefficient>:<Number of independent samples used>)"));
		}

		if (!(trio_sample_info_str.empty())) {

 			cout << "\t - Is sample in corcordant trio (CONC)" << endl;

			if (trio_sample_info_str == "GenomeDKPedigree") {

				all_trio_info = Trio::parseGenomeDKPedigree(vcf_reader.metaData());

			} else {

				all_trio_info = Trio::parsePedigree(vcf_reader.metaData(), trio_sample_info_str);
			}

			output_meta_data.formatDescriptors().emplace("CONC", Attribute::DetailedDescriptor("CONC", Attribute::Number::One, Attribute::Type::String, "Is sample in corcordant trio"));
 		}


        cout << "\n[" << Utils::getLocalTime() << "] Parsing variants ...\n" << endl;

		VcfFileWriter vcf_writer(output_prefix + ".vcf", output_meta_data, true);

		ulong num_variants = 0;
		ulong num_alt_alleles = 0;
		ulong num_trios = 0;

		Variant * cur_var;

        string cur_chromosome = "";

        auto genome_seqs_it = genome_seqs.find(cur_chromosome);
        assert(genome_seqs_it == genome_seqs.end());

		while (vcf_reader.getNextVariant(&cur_var)) {

			num_variants++;

            if (cur_var->chrom() != cur_chromosome) {

                genome_seqs_it = genome_seqs.find(cur_var->chrom());
                cur_chromosome = cur_var->chrom();
            }

			if (!(genome_filename.empty())) {

				assert(genome_seqs_it != genome_seqs.end());
				assert((genome_seqs_it->second->seq().substr(cur_var->pos() - 1, cur_var->ref().seq().size()) == cur_var->ref().seq()) or (cur_var->filters().front() == "UV"));

        		auto homopolymer_info = Auxiliaries::getHomopolymerInfo(cur_var->pos(), genome_seqs_it->second->seq());   
				cur_var->info().setValue<string>("HPL", to_string(homopolymer_info.first) + ":" + homopolymer_info.second);
			}

			if (!(repeat_filename.empty())) {

	            for (uint alt_allele_idx = 0; alt_allele_idx < cur_var->numAlts(); alt_allele_idx++) {

	                num_alt_alleles++;

	                string alt_id = cur_var->chrom() + "_" + to_string(cur_var->pos()) + "_" + to_string(alt_allele_idx);
	                auto repeat_annotations_it = repeat_annotations.find(alt_id);

	                JoiningString alt_rm_label(':');

	                if (repeat_annotations_it != repeat_annotations.end()) {

	                    for (auto & repeat_family : repeat_annotations_it->second) {

	                        JoiningString repeat_info_elements('#');
	                        repeat_info_elements.join(repeat_family.first);
	                        repeat_info_elements.join(to_string(calculateNucleotideCover(repeat_family.second)));

	                        alt_rm_label.join(repeat_info_elements.str());
	                    }

	                    auto allele_attributes = Auxiliaries::alleleAttributes(cur_var->alt(alt_allele_idx), cur_var->ref());
	                    assert((allele_attributes.type == Auxiliaries::Type::Insertion) or (allele_attributes.type == Auxiliaries::Type::Deletion));

	                } else {

	                    alt_rm_label.join(".");
	                }

	                assert(!alt_rm_label.empty());
	                cur_var->alt(alt_allele_idx).info().setValue<string>("RMA", alt_rm_label.str());
	            }
			}

			if (!(indepedent_samples_regex_str.empty())) {

	            auto inbreeding_stats = Stats::calcInbreedingStats(*cur_var, indepedent_samples_regex);

	            if (inbreeding_stats.is_fixed) {

	        		cur_var->info().setValue<string>("IBC", "NA:" + to_string(inbreeding_stats.num_samples));	        	
	            
	            } else {

	            	cur_var->info().setValue<string>("IBC", to_string(abs(inbreeding_stats.inbreeding_coef)) + ":" + to_string(inbreeding_stats.num_samples));
	            }
		    }

			if (!(trio_sample_info_str.empty())) {

				for (auto & trio_info: all_trio_info) {

					num_trios++;

					Trio cur_trio = Trio(*cur_var, trio_info);
					string conc_status = "NA";

					if (!(cur_trio.isFiltered())) {

						if (cur_trio.isInformative()) {

							if (cur_trio.isConcordant()) {

								conc_status = "TRUE";

							} else {

								conc_status = "FALSE";
							}
						}
					}

					cur_var->getSample(trio_info.father).info().setValue<string>("CONC", conc_status);
					cur_var->getSample(trio_info.mother).info().setValue<string>("CONC", conc_status);
					cur_var->getSample(trio_info.child).info().setValue<string>("CONC", conc_status);
				}
			}

            vcf_writer.write(cur_var);
            delete cur_var;

            if ((num_variants % 100000) == 0) {

                std::cout << "[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants" << endl;
            }
		}

		for (auto & genome_seq : genome_seqs) {

            delete genome_seq.second;
        }

		cout << "\n[" << Utils::getLocalTime() << "] Added attributes to " << num_variants << " variants, " << num_alt_alleles << " alternative alleles and " << num_trios << " trios" <<  endl;
		cout << endl;
	}
}
