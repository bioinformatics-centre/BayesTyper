
/*
convertAlleleIdToSequence.cpp - This file is part of BayesTyper (v0.9)


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

#include "vcf++/Utils.hpp"
#include "vcf++/JoiningString.hpp"
#include "vcf++/FastaRecord.hpp"
#include "vcf++/FastaReader.hpp"
#include "vcf++/Auxiliaries.hpp"

#include "convertAlleleIdToSequence.hpp"

const unordered_set<string> relevant_info_att = {"CT", "END", "IMPRECISE", "MEINFO"};

namespace ConvertAlleleIdToSequence {

    unordered_map<string,vector<string> > parseRelevantInfoAttributes(const string & info) {

        unordered_map<string,vector<string> > info_att;

        if (info != ".") {

            auto info_split = Utils::splitString(info, ';');

            for (auto &att: info_split) {

                auto att_split = Utils::splitString(att, '=');
                assert(att_split.size() <= 2);

                if (relevant_info_att.count(att_split.front()) > 0) {

                    if (att_split.size() == 1) {

                        assert(info_att.emplace(att_split.front(), vector<string>()).second);

                    } else {

                        assert(att_split.size() == 2);
                        assert(info_att.emplace(att_split.front(), Utils::splitString(att_split.back(), ',')).second);
                    }
                }
            }
        }

        return info_att;
    }

    uint maxAlleleLength(const string & allele_1, const string & allele_2) {

        auto allele_1_cbit = allele_1.cbegin();
        auto allele_2_cbit = allele_2.cbegin();

        auto allele_1_ceit = allele_1.cend();
        auto allele_2_ceit = allele_2.cend();

        assert(allele_1_cbit != allele_1_ceit);
        assert(allele_2_cbit != allele_2_ceit);

        allele_1_ceit--;
        allele_2_ceit--;

        uint num_trimmed_nt = 0;

        while ((allele_1_cbit != allele_1_ceit) and (allele_2_cbit != allele_2_ceit)) {

            if (*allele_1_ceit != *allele_2_ceit) {

                break;
            }

            num_trimmed_nt++;

            allele_1_ceit--;
            allele_2_ceit--;
        }

        while ((allele_1_cbit != allele_1_ceit) and (allele_2_cbit != allele_2_ceit)) {

            if (*allele_1_cbit != *allele_2_cbit) {

                break;
            }

            num_trimmed_nt++;

            allele_1_cbit++;
            allele_2_cbit++;
        }

        if (*allele_1_cbit == *allele_2_cbit) { 

            num_trimmed_nt++;            
        }

        assert(allele_1.size() >= num_trimmed_nt);
        assert(allele_2.size() >= num_trimmed_nt);

        uint max_length = max(allele_1.size(), allele_2.size()) - num_trimmed_nt;
        assert(max_length > 0);

        return max_length;
    }

    void convertAlleleIdToSequence(const string & vcf_filename, const string & genome_filename, const string & output_prefix, const string & mei_filename, const uint max_allele_length, const bool keep_imprecise_variants, const bool is_delly_output) {

        cout << "[" << Utils::getLocalTime() << "] Running BayesTyperUtils (" << BTU_VERSION << ") convertAlleleIdToSequence ...\n" << endl;

        ifstream in_vcf;
        in_vcf.open(vcf_filename);
        assert(in_vcf.is_open());

        unordered_map<string, FastaRecord*> genome_seqs;
        FastaReader genome_reader(genome_filename);
        FastaRecord * cur_fasta_rec;

        while (genome_reader.getNextRecord(&cur_fasta_rec)) {

            cout << "[" << Utils::getLocalTime() << "] Parsed sequence: " << cur_fasta_rec->id() << endl;
            assert(genome_seqs.emplace(cur_fasta_rec->id(), cur_fasta_rec).second);
        }

        cout << endl;

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

        uint num_vars = 0;
        uint skipped_transversion_variants = 0;
        uint skipped_imprecise_variants = 0;
        uint skipped_seq_match_variants = 0;
        uint skipped_seq_non_canon_variants = 0;

        uint num_alleles = 0;
        uint unsupported_alleles = 0;
        uint max_length_alleles = 0;

        unordered_map<string, int> included_allele_stats;
        unordered_map<string, int> unsupported_allele_stats;
        unordered_map<string, int> max_length_allele_stats;

        regex cnv_regex("^<CN(\\d+)>$");
        regex me_regex("^<INS:ME:(\\S+)>$");
        
        regex tra1_regex("^[\\[|\\]]\\S+:\\d+[\\[|\\]]\\S+$");
        regex tra2_regex("^\\S+[\\[|\\]]\\S+:\\d+[\\[|\\]]$");

        ofstream out_vcf;
        out_vcf.open(output_prefix + ".vcf");

        string in_vcf_line;

        while (getline(in_vcf, in_vcf_line)) {

            if (in_vcf_line.front() == '#') {

                if (!(in_vcf_line.substr(0,6) == "##INFO")) {

                    if (in_vcf_line.substr(0,8) == "##FORMAT") {

                        if (in_vcf_line.find("<ID=GT,") != string::npos) {

                            out_vcf << in_vcf_line << "\n";
                        }

                    } else {

                        out_vcf << in_vcf_line << "\n";
                    }
                } 

            } else {

                num_vars++;

                vector<string> in_vcf_line_split = Utils::splitString(in_vcf_line, '\t');
                assert(in_vcf_line_split.size() >= 8);

                auto chr_seq = genome_seqs.find(in_vcf_line_split.at(0));
                assert(chr_seq != genome_seqs.end());

                assert(!(in_vcf_line_split.at(3).empty()));

                const int var_pos = stoi(in_vcf_line_split.at(1));
                auto relevant_info_att = parseRelevantInfoAttributes(in_vcf_line_split.at(7));

                if ((in_vcf_line_split.at(4) == "<TRA>") or (regex_match(in_vcf_line_split.at(4), tra1_regex)) or (regex_match(in_vcf_line_split.at(4), tra2_regex))) {

                    skipped_transversion_variants++;
                    continue;
                } 

                if (!keep_imprecise_variants) {

                    if (relevant_info_att.count("IMPRECISE") > 0) {

                        skipped_imprecise_variants++;
                        continue;
                    }
                }

                auto relevant_info_att_ct_it = relevant_info_att.find("CT");

                if (is_delly_output) {

                    assert(in_vcf_line_split.at(3) == "N");

                    assert(relevant_info_att_ct_it != relevant_info_att.end());
                    assert(relevant_info_att_ct_it->second.size() == 1);
                }

                if (in_vcf_line_split.at(3) == "N") {

                    in_vcf_line_split.at(3) = chr_seq->second->seq().at(var_pos - 1);
                }

                const string anchor_nt = in_vcf_line_split.at(3).substr(0,1);
                string var_seq = "";

                if (chr_seq->second->seq().substr(var_pos - 1, in_vcf_line_split.at(3).size()) != in_vcf_line_split.at(3)) {
                
                    skipped_seq_match_variants++;
                    continue;
                }

                auto relevant_info_att_end_it = relevant_info_att.find("END");

                if (relevant_info_att_end_it != relevant_info_att.end()) {

                    assert(relevant_info_att_end_it->second.size() == 1);
                    int var_end = stoi(relevant_info_att_end_it->second.front());

                    if (is_delly_output) {

                        if (relevant_info_att_ct_it->second.front().back() == '5') {

                            var_end--;
                        }
                    }

                    assert(var_end > 0);
                    assert(var_pos <= var_end);

                    int ref_len = var_end - var_pos;
                    var_seq = chr_seq->second->seq().substr(var_pos, ref_len);
                }

                if ((in_vcf_line_split.at(3).find_first_not_of("ACGTN") != string::npos) or (var_seq.find_first_not_of("ACGTN") != string::npos)) {

                    skipped_seq_non_canon_variants++;
                    continue;
                }

                auto alt_split = Utils::splitString(in_vcf_line_split.at(4), ',');
                JoiningString alt_elem(',');

                for (uint alt_idx = 0; alt_idx < alt_split.size(); alt_idx++) {

                    num_alleles++;

                    smatch cnv_match_res;
                    smatch me_match_res;
                    smatch ins_match_res;

                    string alt_seq = ""; 
                    string allele_type = "";

                    if (alt_split.at(alt_idx).front() != '<') {

                        assert(!is_delly_output);
                        assert(alt_split.at(alt_idx).back() != '<');

                        if (alt_split.at(alt_idx).find_first_not_of("ACGTN") != string::npos) {

                            allele_type = "non-canonical base(s)";

                        } else {

                            alt_seq += alt_split.at(alt_idx);

                            if (!var_seq.empty()) { // Needed because sequence alts not guaranteed to have END coordinate

                                assert(var_seq.size() + 1 == in_vcf_line_split.at(3).size());
                                assert(var_seq == chr_seq->second->seq().substr(var_pos, in_vcf_line_split.at(3).size() - 1));
                            }                       
    
                            allele_type = "sequence";
                        }                        

                    } else if (regex_match(alt_split.at(alt_idx), cnv_match_res, cnv_regex)) {

                        assert(!var_seq.empty());
                        assert(in_vcf_line_split.at(3).size() == 1);
                        assert(!is_delly_output);

                        in_vcf_line_split.at(3) = anchor_nt + var_seq;

                        assert(cnv_match_res.size() == 2);
                        int copy_number = stoi(cnv_match_res[1].str());
                        assert(copy_number >= 0);

                        if (copy_number == 0) { 

                            alt_seq += anchor_nt;

                        } else {

                            assert(copy_number > 1);
                            alt_seq += anchor_nt + var_seq;

                            for (int i = 1; i < copy_number; ++i) {

                                alt_seq += var_seq;
                            }
                        }

                        allele_type = alt_split.at(alt_idx);

                    } else if ((alt_split.at(alt_idx) == "<DUP>") or (alt_split.at(alt_idx) == "<DUP:TANDEM>")) {

                        assert(!var_seq.empty());
                        assert(in_vcf_line_split.at(3).size() == 1);
                        assert(alt_split.size() == 1);

                        if (is_delly_output) {

                            assert(relevant_info_att_ct_it != relevant_info_att.end());
                            assert(relevant_info_att_ct_it->second.front() == "5to3");
                            
                            in_vcf_line_split.at(3) = anchor_nt + var_seq;
                            alt_seq += anchor_nt + var_seq + anchor_nt + var_seq;
                        
                        } else {

                            in_vcf_line_split.at(3) = anchor_nt + var_seq;
                            alt_seq += anchor_nt + var_seq + var_seq;
                        }

                        allele_type = alt_split.at(alt_idx);

                    } else if (alt_split.at(alt_idx) == "<DEL>") {

                        assert(!var_seq.empty());
                        assert(in_vcf_line_split.at(3).size() == 1);
                        assert(alt_split.size() == 1);

                        if (is_delly_output) {

                            assert(relevant_info_att_ct_it != relevant_info_att.end());
                            assert(relevant_info_att_ct_it->second.front() == "3to5");
                        }
                        
                        in_vcf_line_split.at(3) = anchor_nt + var_seq;
                        alt_seq += anchor_nt;

                        allele_type = alt_split.at(alt_idx);

                    } else if (alt_split.at(alt_idx) == "<INV>") {

                        assert(!var_seq.empty());
                        assert(in_vcf_line_split.at(3).size() == 1);
                        assert(alt_split.size() == 1);

                        if (is_delly_output) {

                            assert(relevant_info_att_ct_it != relevant_info_att.end());

                            if (relevant_info_att_ct_it->second.front() == "3to3") {

                                in_vcf_line_split.at(3) = anchor_nt + var_seq;
                                alt_seq += anchor_nt + Auxiliaries::reverseComplementSequence(var_seq);                                
                            
                            } else {
                                
                                assert(relevant_info_att_ct_it->second.front() == "5to5");

                                in_vcf_line_split.at(3) = anchor_nt + var_seq;
                                alt_seq += Auxiliaries::reverseComplementSequence(anchor_nt + var_seq);
                            }

                        } else {

                            in_vcf_line_split.at(3) = anchor_nt + var_seq;
                            alt_seq += anchor_nt + Auxiliaries::reverseComplementSequence(var_seq);
                        }

                        allele_type = alt_split.at(alt_idx);

                    } else if (regex_match(alt_split.at(alt_idx), me_match_res, me_regex)) {

                        assert(in_vcf_line_split.at(3).size() == 1);
                        assert(alt_split.size() == 1);
                        assert(!is_delly_output);

                        assert(me_match_res.size() == 2);
                        auto me_type = me_match_res[1].str();

                        auto me_seqs_it = mei_seqs.find(me_type);

                        auto relevant_info_att_meinfo_it = relevant_info_att.find("MEINFO");

                        if ((me_seqs_it != mei_seqs.end()) and (relevant_info_att_meinfo_it != relevant_info_att.end())) {
                            
                            assert(relevant_info_att_meinfo_it->second.size() == 4);

                            transform(me_type.begin(), me_type.end(), me_type.begin(), ::toupper);
                            transform(relevant_info_att_meinfo_it->second.front().begin(), relevant_info_att_meinfo_it->second.front().end(), relevant_info_att_meinfo_it->second.front().begin(), ::toupper);

                            assert(relevant_info_att_meinfo_it->second.front().find(me_type) != string::npos);

                            FastaRecord * me_type_seq = nullptr;

                            if (relevant_info_att_meinfo_it->second.back() == "+") {

                                me_type_seq = me_seqs_it->second.first;

                            } else {

                                assert(relevant_info_att_meinfo_it->second.back() == "-");
                                me_type_seq = me_seqs_it->second.second;
                            }

                            assert(me_type_seq);

                            if (stoi(relevant_info_att_meinfo_it->second.at(1)) > 0) {

                                assert(static_cast<uint>(stoi(relevant_info_att_meinfo_it->second.at(2))) <= me_type_seq->seq().size());
                                
                                alt_seq += anchor_nt + string(me_type_seq->seq().begin() + stoi(relevant_info_att_meinfo_it->second.at(1)) - 1, me_type_seq->seq().begin() + stoi(relevant_info_att_meinfo_it->second.at(2)));

                                allele_type = alt_split.at(alt_idx);

                            } else {

                                assert(stoi(relevant_info_att_meinfo_it->second.at(1)) == 0);
                                assert(stoi(relevant_info_att_meinfo_it->second.at(2)) == 0);
                            }
                        }
                    }

                    if (alt_seq.empty()) {

                        if (allele_type.empty()) {

                            assert(alt_split.at(alt_idx).front() == '<');
                            assert(alt_split.at(alt_idx).back() == '>');

                            allele_type = alt_split.at(alt_idx);
                        }
                        
                        assert(alt_split.size() == 1);
                        unsupported_alleles++;

                        auto unsupported_allele_stats_it = unsupported_allele_stats.emplace(allele_type, 0);
                        unsupported_allele_stats_it.first->second++;
                    
                    } else if (maxAlleleLength(in_vcf_line_split.at(3), alt_seq) > max_allele_length) {

                        assert(!(allele_type.empty()));

                        max_length_alleles++;

                        auto max_length_allele_stats_it = max_length_allele_stats.emplace(allele_type, 0);
                        max_length_allele_stats_it.first->second++;                        

                    } else {

                        assert(!(allele_type.empty()));

                        alt_elem.join(alt_seq);

                        auto included_allele_stats_it = included_allele_stats.emplace(allele_type, 0);
                        included_allele_stats_it.first->second++;
                    }
                }

                if (!(alt_elem.empty())) {

                    in_vcf_line_split.at(4) = alt_elem.str();
                    in_vcf_line_split.at(7) = ".";

                    if (in_vcf_line_split.size() > 8) {

                        auto format_split = Utils::splitString(in_vcf_line_split.at(8), ':');
                        assert(format_split.front() == "GT");

                        in_vcf_line_split.at(8) = "GT";

                        for (uint i = 9; i < in_vcf_line_split.size(); i++) {

                            in_vcf_line_split.at(i) = Utils::splitString(in_vcf_line_split.at(i), ':').front();
                        }
                    }
                    
                    JoiningString out_vcf_line('\t');
                    out_vcf_line.join(in_vcf_line_split);

                    out_vcf << out_vcf_line.str() + "\n";
                }

                if (num_vars % 100000 == 0) {

                    cout << "[" << Utils::getLocalTime() << "] Parsed " << num_vars << " variant(s) ..." << endl;
                }
            }
        }

        in_vcf.close();
        out_vcf.close();

        for (auto & genome_seq : genome_seqs) {

            delete genome_seq.second;
        }

        for (auto & mei_seq : mei_seqs) {

            delete mei_seq.second.first;
            delete mei_seq.second.second;
        }

        cout << "\n[" << Utils::getLocalTime() << "] Completed parsing of " << num_vars << " variant(s)" << endl;

        cout << "\n\t- Skipped " << skipped_transversion_variants << " transversion variants(s)." << endl;

        if (!keep_imprecise_variants) {

            cout << "\t- Skipped " << skipped_imprecise_variants << " imprecise variants(s)." << endl;            
        }

        cout << "\n\t- Skipped " << skipped_seq_match_variants << " variants(s) due to reference - vcf sequence mismatch." << endl;
        cout << "\t- Skipped " << skipped_seq_non_canon_variants << " variants(s) due to reference sequence having non-canonical base(s)." << endl;

        cout << "\n[" << Utils::getLocalTime() << "] Parsed " << num_alleles << " allele(s)\n" << endl;

        for (auto & key_val: included_allele_stats) {

            cout << "\t- Included " << key_val.second << " " << key_val.first << " allele(s) " << endl;
        }

        cout << "\n\t- Skipped " << unsupported_alleles << " unsupported alleles:\n" << endl;

        uint unsupported_allele_sum = 0;
 
        for (auto & key_val : unsupported_allele_stats) {

            unsupported_allele_sum += key_val.second;
            cout << "\t\t- " << key_val.second << " " << key_val.first << " allele(s)" << endl;
        }

        assert(unsupported_allele_sum == unsupported_alleles);

        cout << "\n\t- Skipped " << max_length_alleles << " allele(s) larger than " << max_allele_length << " bases:\n" << endl;

        uint max_length_allele_sum = 0;
 
        for (auto & key_val : max_length_allele_stats) {

            max_length_allele_sum += key_val.second;
            cout << "\t\t- " << key_val.second << " " << key_val.first << " allele(s)" << endl;
        }

        assert(max_length_allele_sum == max_length_alleles);

        cout << endl;
    }
}

