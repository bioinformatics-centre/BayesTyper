
/*
addRepeatAnnotation.cpp - This file is part of BayesTyper (v0.9)


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


#include <fstream>
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
#include <algorithm>

#include "vcf++/Utils.hpp"
#include "vcf++/VcfFile.hpp"
#include "vcf++/Variant.hpp"
#include "vcf++/JoiningString.hpp"
#include "vcf++/Auxiliaries.hpp"

#include "addRepeatAnnotation.hpp"

using namespace std;

namespace AddRepeatAnnotation {

    struct SortCover {
        
        bool operator() (const pair<uint, uint> & lhs, const pair<uint, uint> & rhs) {

            return lhs.first < rhs.first;
        }
    };

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

            } 

            cover_it++;
        }

        return nt_cover;
    }

    void addRepeatAnnotation(const string & in_vcf_filename, const string & repeat_masker_filename, const string & output_prefix) {

        cout << "[" << Utils::getLocalTime() << "] Running BayesTyperUtils (" << BTU_VERSION << ") addRepeatAnnotation ...\n" << endl;

        cout << "[" << Utils::getLocalTime() << "] Parsing RepeatMasker annotation ..." << endl;

        ifstream rm_file(repeat_masker_filename);

        unordered_map<string, map<string, vector<pair<uint, uint> > > > rm_allele_annotations;

        string line;
        uint line_count = 0;

        while (getline(rm_file, line)) {

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

                auto rm_allele_annotations_it = rm_allele_annotations.emplace(line_split.at(4), map<string, vector<pair<uint, uint> > >());
                auto rm_allele_annotations_family_it = rm_allele_annotations_it.first->second.emplace(name, vector<pair<uint, uint> >());

                rm_allele_annotations_family_it.first->second.emplace_back(query_start, query_end);
            }
        }

        rm_file.close();

        map<uint,uint> label_count_freqs;

        for (auto & rm_allele_annotation : rm_allele_annotations) {

            uint num_labels = rm_allele_annotation.second.size();

            auto label_count_freq_emplace_res = label_count_freqs.emplace(num_labels, 1);

            if (!label_count_freq_emplace_res.second) {

                label_count_freq_emplace_res.first->second++;
            }
        }

        cout << "[" << Utils::getLocalTime() << "] Completed parsing of RepeatMasker annotation \n" << endl;
        cout << "[" << Utils::getLocalTime() << "] Label number distribution:\n" << endl;

        for (auto & label_count_freq : label_count_freqs) {

            cout << "\t- " << label_count_freq.first << ": " << label_count_freq.second << endl;
        }

        cout << "\n[" << Utils::getLocalTime() << "] Adding RepeatMasker labels to vcf ...\n" << endl;

        GenotypedVcfFileReader in_vcf(in_vcf_filename, true);

        VcfMetaData out_meta_data = in_vcf.metaData();
        assert(out_meta_data.infoDescriptors().emplace("RMA", Attribute::DetailedDescriptor("RMA", Attribute::Number::A, Attribute::Type::String, "RepeatMasker annotations (<family#nucleotide_cover>:...)")).second);

        VcfFileWriter out_vcf(output_prefix + ".vcf", out_meta_data, true);

        Variant * cur_var;
        uint num_variants = 0;
        uint num_allelles = 0;
        uint num_allelles_labelled = 0;

        while (in_vcf.getNextVariant(&cur_var)) {

            num_variants++;

            for (uint alt_idx = 0; alt_idx < cur_var->numAlts(); alt_idx++) {

                num_allelles++;

                string alt_id = cur_var->chrom() + "_" + to_string(cur_var->pos()) + "_" + to_string(alt_idx);
                auto rm_allele_annotations_find_res = rm_allele_annotations.find(alt_id);

                JoiningString alt_rm_label(':');

                if (rm_allele_annotations_find_res != rm_allele_annotations.end()) {

                    num_allelles_labelled++;

                    for (auto &repeat_family : rm_allele_annotations_find_res->second) {

                        JoiningString repeat_info_elements('#');
                        repeat_info_elements.join(repeat_family.first);
                        repeat_info_elements.join(to_string(calculateNucleotideCover(repeat_family.second)));

                        alt_rm_label.join(repeat_info_elements.str());
                    }

                    auto allele_attributes = Auxiliaries::alleleAttributes(cur_var->alt(alt_idx), cur_var->ref());
                    assert((allele_attributes.type == Auxiliaries::Type::Insertion) or (allele_attributes.type == Auxiliaries::Type::Deletion));

                } else {

                    alt_rm_label.join(".");
                }

                assert(!alt_rm_label.empty());
                assert(cur_var->alt(alt_idx).info().setValue("RMA", alt_rm_label.str()));
            }

            out_vcf.write(cur_var);

            if (num_variants % 100000 == 0) {

                cout << "[" << Utils::getLocalTime() << "] Labelled " << num_variants << " variant(s)" << endl;
            }

            delete cur_var;
        }

        cout << "\n[" << Utils::getLocalTime() << "] Completed BayesTyperUtils addRepeatMaskerAnnotation" << endl;
        cout << "[" << Utils::getLocalTime() << "] " << num_variants << " variant(s) were parsed in total" << endl;
        cout << "[" << Utils::getLocalTime() << "] " << num_allelles << " allele(s) were parsed in total" << endl;
        cout << "[" << Utils::getLocalTime() << "] " << num_allelles_labelled << " allele(s) were labelled\n" << endl;
    }
}
