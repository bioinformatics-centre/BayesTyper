
/*
convertAlleleId.cpp - This file is part of BayesTyper (v0.9)


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
#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>
#include <regex>
#include <unordered_set>
#include <algorithm>

#include "Utils.hpp"
#include "JoiningString.hpp"
#include "FastaRecord.hpp"
#include "FastaReader.hpp"
#include "Auxiliaries.hpp"

#include "convertAlleleId.hpp"


namespace ConvertAlleleId {

    void convertAlleleId(const string & vcf_filename, const string & genome_filename, const string & output_prefix, const string & mei_filename, const bool keep_imprecise_variants) {

        cout << "[" << Utils::getLocalTime() << "] Running BayesTyperTools (" << BT_VERSION << ") convertAlleleId ...\n" << endl;

        unordered_map<string, FastaRecord*> genome_seqs;
        FastaReader genome_reader(genome_filename);
        FastaRecord * cur_fasta_rec;

        while (genome_reader.getNextRecord(&cur_fasta_rec)) {

            assert(genome_seqs.emplace(cur_fasta_rec->id(), cur_fasta_rec).second);
            genome_seqs.at(cur_fasta_rec->id())->convertToUppercase();
        }

        cout << "[" << Utils::getLocalTime() << "] Parsed " << genome_seqs.size() << " chromosome(s)" << endl;

        unordered_map<string, pair<FastaRecord*, FastaRecord*> > mei_seqs;

        if (!(mei_filename.empty())) {

            FastaReader mei_reader(mei_filename);

            while (mei_reader.getNextRecord(&cur_fasta_rec)) {

                FastaRecord * cur_fasta_rec_rv = new FastaRecord(cur_fasta_rec->id(), Auxiliaries::reverseComplementSequence(cur_fasta_rec->seq()));

                cur_fasta_rec->convertToUppercase();
                cur_fasta_rec_rv->convertToUppercase();

                assert(mei_seqs.emplace(cur_fasta_rec->id(), make_pair(cur_fasta_rec, cur_fasta_rec_rv)).second);
            }
        } 

        cout << "[" << Utils::getLocalTime() << "] Parsed " << mei_seqs.size() << " mobile element insertion sequence(s)\n" << endl;

        uint num_variants = 0;
        uint skipped_imprecise_variants = 0;
        uint skipped_seq_match_variants = 0;

        uint num_alt_alleles = 0;
        uint num_unsupported_alt_alleles = 0;

        unordered_map<string, int> included_alt_allele_stats;
        unordered_map<string, int> unsupported_alt_allele_stats;

        regex cnv_regex("^<CN(\\d+)>$");
        regex me_regex("^<INS:ME:(\\S+)>$");
        
        regex tra1_regex("^[\\[|\\]]\\S+:\\d+[\\[|\\]]\\S+$");
        regex tra2_regex("^\\S+[\\[|\\]]\\S+:\\d+[\\[|\\]]$");

        GenotypedVcfFileReader vcf_reader(vcf_filename, false);

        auto output_meta_data = vcf_reader.metaData();
        output_meta_data.formatDescriptors().clear();

        VcfFileWriter vcf_writer(output_prefix + ".vcf", output_meta_data, false);
        
        Variant * cur_var;
        string cur_chromosome = "";

        auto genome_seqs_it = genome_seqs.find(cur_chromosome);
        assert(genome_seqs_it == genome_seqs.end());

        while (vcf_reader.getNextVariant(&cur_var)) {

            num_variants++;

            if (cur_var->chrom() != cur_chromosome) {

                genome_seqs_it = genome_seqs.find(cur_var->chrom());
                assert(genome_seqs_it != genome_seqs.end());

                cur_chromosome = cur_var->chrom();
            }

            if (!keep_imprecise_variants) {

                if (cur_var->info().getValue<string>("IMPRECISE").second) {

                    skipped_imprecise_variants++;
                    delete cur_var;

                    continue;
                }
            }

            if (cur_var->ref().seq() == "N") {

                cur_var->ref().seq() = genome_seqs_it->second->seq().at(cur_var->pos() - 1);
            }

            if (genome_seqs_it->second->seq().substr(cur_var->pos() - 1, cur_var->ref().seq().size()) != cur_var->ref().seq()) {
            
                skipped_seq_match_variants++;
                delete cur_var;

                continue;
            }

            const string anchor_nt = cur_var->ref().seq().substr(0, 1);
            string end_variant_seq = "";

            auto svlen_value = cur_var->info().getValue<int>("SVLEN");

            if (svlen_value.second) {

                auto end_value = cur_var->info().getValue<int>("END");

                if (end_value.second) {

                    assert(end_value.first > 0);
                    assert(static_cast<uint>(end_value.first) == (cur_var->pos() + abs(svlen_value.first)));

                } else {

                    assert(cur_var->info().setValue<int>("END", cur_var->pos() + abs(svlen_value.first)));
                }
            }

            auto end_value = cur_var->info().getValue<int>("END");

            if (end_value.second) {

                assert(end_value.first > 0);

                if (cur_var->pos() <= static_cast<uint>(end_value.first)) {

                    end_variant_seq = genome_seqs_it->second->seq().substr(cur_var->pos(), end_value.first - cur_var->pos());
                }
            } 

            vector<uint> excluded_alt_indices;

            bool variant_has_sequence = false;
            bool variant_has_id = false;

            const uint cur_num_alts = cur_var->numAlts();

            for (uint alt_allele_idx = 0; alt_allele_idx < cur_num_alts; alt_allele_idx++) {

                num_alt_alleles++;

                smatch cnv_match_res;
                smatch me_match_res;
                smatch ins_match_res;

                bool is_excluded = false;
                string allele_type = "";

                if (cur_var->alt(alt_allele_idx).isID()) {

                    assert(!variant_has_sequence);
                    variant_has_id = true;

                    allele_type = cur_var->alt(alt_allele_idx).seq();

                    if ((cur_var->alt(alt_allele_idx).seq() == "<TRA>") or (regex_match(cur_var->alt(alt_allele_idx).seq(), tra1_regex)) or (regex_match(cur_var->alt(alt_allele_idx).seq(), tra2_regex))) {

                        allele_type = "transversion";

                        is_excluded = true;
                        excluded_alt_indices.push_back(alt_allele_idx);

                    } else if (cur_var->alt(alt_allele_idx).seq() == "<CNV>") {

                        assert(cur_var->numAlts() == 1);

                        auto gscndist_value = cur_var->info().getValue<string>("GSCNDIST");
                        assert(gscndist_value.second);

                        auto gscndist_value_split = Utils::splitString(gscndist_value.first, ','); 
                        assert(gscndist_value_split.size() > 1);
                        
                        vector<uint> cnv_multiplicities;
                        cnv_multiplicities.reserve(gscndist_value_split.size());

                        for (uint i = 0; i < gscndist_value_split.size(); i++) {

                            if (i == 1) {

                                continue;
                            }

                            if (stoi(gscndist_value_split.at(i)) > 0) {

                                cnv_multiplicities.push_back(i);
                            }
                        } 

                        if (cnv_multiplicities.empty() or (end_variant_seq.empty())) {

                            is_excluded = true;
                            excluded_alt_indices.push_back(alt_allele_idx);
                        
                        } else {

                            cur_var->ref().seq() = anchor_nt + end_variant_seq;
                            
                            cur_var->removeAlts({0});
                            assert(cur_var->numAlts() == 0);
                                
                            for (auto & cnv_multiplicity: cnv_multiplicities) {

                                if (cnv_multiplicity == 0) { 

                                    assert(cur_var->addAlt(Allele(anchor_nt)).second);

                                } else {

                                    assert(cnv_multiplicity > 1);
                                    string cnv_seq = anchor_nt + end_variant_seq;

                                    for (uint i = 1; i < cnv_multiplicity; i++) {

                                        cnv_seq += end_variant_seq;
                                    }

                                    assert(cur_var->addAlt(Allele(cnv_seq)).second);                                    
                                }
                            }
                        } 

                    } else if (regex_match(cur_var->alt(alt_allele_idx).seq(), cnv_match_res, cnv_regex)) {

                        if (end_variant_seq.empty()) {

                            is_excluded = true;
                            excluded_alt_indices.push_back(alt_allele_idx);

                        } else {

                            cur_var->ref().seq() = anchor_nt + end_variant_seq;

                            assert(cnv_match_res.size() == 2);
                            int cnv_multiplicity = stoi(cnv_match_res[1].str());
                            assert(cnv_multiplicity >= 0);

                            if (cnv_multiplicity == 0) { 

                                cur_var->alt(alt_allele_idx).seq() = anchor_nt;

                            } else if (cnv_multiplicity > 1) {

                                string cnv_seq = anchor_nt + end_variant_seq;

                                for (int i = 1; i < cnv_multiplicity; ++i) {

                                    cnv_seq += end_variant_seq;
                                }

                                cur_var->alt(alt_allele_idx).seq() = cnv_seq;                                
                            
                            } else {

                                is_excluded = true;
                                excluded_alt_indices.push_back(alt_allele_idx);                                
                            }
                        }

                    } else if ((cur_var->alt(alt_allele_idx).seq() == "<DUP>") or (cur_var->alt(alt_allele_idx).seq() == "<DUP:")) {

                        if (end_variant_seq.empty()) {

                            is_excluded = true;
                            excluded_alt_indices.push_back(alt_allele_idx);

                        } else {

                            cur_var->ref().seq() = anchor_nt + end_variant_seq;
                            cur_var->alt(alt_allele_idx).seq() = anchor_nt + end_variant_seq + end_variant_seq;                                
                        }

                    } else if ((cur_var->alt(alt_allele_idx).seq() == "<DEL>") or (cur_var->alt(alt_allele_idx).seq().substr(0,5) == "<DEL:")) {

                        if (end_variant_seq.empty()) {

                            is_excluded = true;
                            excluded_alt_indices.push_back(alt_allele_idx);

                        } else {

                            cur_var->ref().seq() = anchor_nt + end_variant_seq;
                            cur_var->alt(alt_allele_idx).seq() = anchor_nt;                                
                        }

                    } else if (cur_var->alt(alt_allele_idx).seq() == "<INV>") {

                        if (end_variant_seq.empty() or (end_variant_seq.find_first_not_of("ACGTN") != string::npos)) {

                            is_excluded = true;
                            excluded_alt_indices.push_back(alt_allele_idx);

                        } else {

                            cur_var->ref().seq() = anchor_nt + end_variant_seq;
                            cur_var->alt(alt_allele_idx).seq() = anchor_nt + Auxiliaries::reverseComplementSequence(end_variant_seq);                                
                        }

                    } else if (regex_match(cur_var->alt(alt_allele_idx).seq(), me_match_res, me_regex)) {

                        assert(cur_var->numAlts() == 1);

                        assert(cur_var->ref().seq() == anchor_nt);

                        assert(me_match_res.size() == 2);
                        auto me_type = me_match_res[1].str();

                        auto me_seqs_it = mei_seqs.find(me_type);
                        auto meinfo_value = cur_var->info().getValue<string>("MEINFO");

                        if ((me_seqs_it != mei_seqs.end()) and meinfo_value.second) {

                            auto meinfo_value_split = Utils::splitString(meinfo_value.first, ','); 
                            assert(meinfo_value_split.size() == 4);

                            transform(me_type.begin(), me_type.end(), me_type.begin(), ::toupper);
                            transform(meinfo_value_split.front().begin(), meinfo_value_split.front().end(), meinfo_value_split.front().begin(), ::toupper);

                            assert(meinfo_value_split.front().find(me_type) != string::npos);

                            FastaRecord * me_type_seq = nullptr;

                            if (meinfo_value_split.back() == "+") {

                                me_type_seq = me_seqs_it->second.first;

                            } else {

                                assert(meinfo_value_split.back() == "-");
                                me_type_seq = me_seqs_it->second.second;
                            }

                            assert(me_type_seq);

                            if (stoi(meinfo_value_split.at(1)) > 0) {

                                assert(static_cast<uint>(stoi(meinfo_value_split.at(2))) <= me_type_seq->seq().size());
                                
                                cur_var->alt(alt_allele_idx).seq() = anchor_nt + string(me_type_seq->seq().begin() + stoi(meinfo_value_split.at(1)) - 1, me_type_seq->seq().begin() + stoi(meinfo_value_split.at(2)));

                            } else {

                                assert(stoi(meinfo_value_split.at(1)) == 0);
                                assert(stoi(meinfo_value_split.at(2)) == 0);

                                is_excluded = true;
                                excluded_alt_indices.push_back(alt_allele_idx);
                            }
                        
                        } else {

                            is_excluded = true;
                            excluded_alt_indices.push_back(alt_allele_idx);
                        }
                    
                    } else {

                        is_excluded = true;
                        excluded_alt_indices.push_back(alt_allele_idx);                        
                    }

                } else {

                    assert(!variant_has_id);
                    variant_has_sequence = true;

                    if (cur_var->alt(alt_allele_idx).isMissing()) {
                        
                        allele_type = "missing";

                    } else {

                        allele_type = "sequence";                        
                    }
                }

                assert(!(allele_type.empty()));

                if (is_excluded) {
                    
                    num_unsupported_alt_alleles++;

                    auto unsupported_alt_allele_stats_it = unsupported_alt_allele_stats.emplace(allele_type, 0);
                    unsupported_alt_allele_stats_it.first->second++;             

                } else {

                    auto included_alt_allele_stats_it = included_alt_allele_stats.emplace(allele_type, 0);
                    included_alt_allele_stats_it.first->second++;
                }
            }

            assert(cur_var->numAlts() > 0);
            assert(cur_var->numAlts() >= excluded_alt_indices.size());

            if (!(excluded_alt_indices.empty())) {

                cur_var->removeAlts(excluded_alt_indices);
            }

            if (cur_var->numAlts() > 0) {

                vcf_writer.write(cur_var);
            }

            delete cur_var;

            if (num_variants % 100000 == 0) {

                cout << "[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variant(s) ..." << endl;
            }
        }

        for (auto & genome_seq : genome_seqs) {

            delete genome_seq.second;
        }

        for (auto & mei_seq : mei_seqs) {

            delete mei_seq.second.first;
            delete mei_seq.second.second;
        }

        cout << "\n[" << Utils::getLocalTime() << "] Completed parsing of " << num_variants << " variant(s)" << endl;


        cout << "\n\t- Skipped " << skipped_imprecise_variants << " imprecise variants(s)." << endl;            
        cout << "\t- Skipped " << skipped_seq_match_variants << " variants(s) due to reference - vcf sequence mismatch." << endl;

        cout << "\n[" << Utils::getLocalTime() << "] Parsed " << num_alt_alleles << " alternative allele(s)\n" << endl;

        for (auto & included_alt_allele: included_alt_allele_stats) {

            cout << "\t- Included " << included_alt_allele.second << " " << included_alt_allele.first << " alternative allele(s) " << endl;
        }

        cout << "\n\t- Skipped " << num_unsupported_alt_alleles << " unsupported allele(s):\n" << endl;

        uint unsupported_alt_allele_sum = 0;
 
        for (auto & unsupported_alt_allele : unsupported_alt_allele_stats) {

            unsupported_alt_allele_sum += unsupported_alt_allele.second;
            cout << "\t\t- " << unsupported_alt_allele.second << " " << unsupported_alt_allele.first << " alternative allele(s)" << endl;
        }

        assert(unsupported_alt_allele_sum == num_unsupported_alt_alleles);

        cout << endl;
    }
}

