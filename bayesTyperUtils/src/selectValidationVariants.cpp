
/*
selectValidationVariants.cpp - This file is part of BayesTyper (v0.9)


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
#include <random>

#include "vcf++/VcfFile.hpp"
#include "vcf++/Variant.hpp"
#include "vcf++/FastaRecord.hpp"
#include "vcf++/FastaReader.hpp"
#include "vcf++/JoiningString.hpp"
#include "vcf++/Stats.hpp"
#include "vcf++/Auxiliaries.hpp"
#include "vcf++/Contig.hpp"

#include "selectValidationVariants.hpp"

namespace SelectValidationVariants {

    void selectValidationVariants(const string & in_vcf_filename, const vector<string> & val_trio_ids, const string & genome_fasta_filename, const string & output_prefix, int min_flank_length, int max_fragment_length, float min_gpp, float min_acp) {

        cout << "[" << Utils::getLocalTime() << "] Running BayesTyperUtils (" << BTU_VERSION << ") selectValidationVariants ...\n" << endl;

        vector<string> allele_types = {"SNP", "Insertion", "Deletion", "Inversion", "Complex"};

        unordered_map<string, unique_ptr<ofstream> > output_files;
        unordered_map<string, unique_ptr<ofstream> > output_files_denovo;

        for (auto & allele_type : allele_types) {

            output_files.emplace(allele_type, unique_ptr<ofstream>(new ofstream(output_prefix + "_" + allele_type + ".txt")));
            output_files_denovo.emplace(allele_type, unique_ptr<ofstream>(new ofstream(output_prefix + "_" + allele_type + "_denovo.txt")));
        }

        string header("chr\tpos\tmax_all_length\tmax_all_trim_length\tmax_het_gt_post_sample_id\tmax_get_gt_post\tupstream_flank_start\tupstream_flank_end\tdownstream_flank_start\tdownstream_flank_end\tfastg_seq\talt_ac\tvar_an");

        for (auto & out_file : output_files) {

            *out_file.second << header << endl;
        }

        for (auto & out_file_denovo : output_files_denovo) {

            *out_file_denovo.second << header << endl;
        }

        cout << "[" << Utils::getLocalTime() << "] Parsing reference genome fasta ...\n" << endl;

        unordered_map<string,FastaRecord*> genome_seqs;
        FastaReader genome_reader(genome_fasta_filename);
        FastaRecord * cur_fasta_rec;

        while (genome_reader.getNextRecord(&cur_fasta_rec)) {

            assert(genome_seqs.emplace(cur_fasta_rec->id(), cur_fasta_rec).second);
        }

        GenotypedVcfFileReader in_vcf(in_vcf_filename, true);

		in_vcf.metaData().infoDescriptors().erase("ACO");
        in_vcf.metaData().infoDescriptors().erase("AsmVar_ASQR");

        assert(in_vcf.metaData().formatDescriptors().erase("MAP"));
        assert(in_vcf.metaData().formatDescriptors().erase("NK"));
        assert(in_vcf.metaData().formatDescriptors().erase("NOK"));
        assert(in_vcf.metaData().formatDescriptors().erase("NUK"));

        auto sample_ids = in_vcf.metaData().sampleIds();

        if (!val_trio_ids.empty()) {

            uint num_rm_sample_ids = 0;

            for (auto & sample_id : sample_ids) {

                if (find(val_trio_ids.begin(), val_trio_ids.end(), Utils::splitString(sample_id, '-').front()) == val_trio_ids.end()) {

                    in_vcf.metaData().rmSample(sample_id);
                    num_rm_sample_ids++;
                }
            }

            cout << "[" << Utils::getLocalTime() << "] Removed " << num_rm_sample_ids << " sample id(s)\n" << endl;
            sample_ids = in_vcf.metaData().sampleIds();
        }

        int upstream_var_end = 0;
        int downstream_var_start = 0;

        Variant * cur_eval_var = nullptr;
        Variant * cur_downstream_var;

        string cur_eval_chr;
        string cur_downstream_var_chrom;

        bool downstream_var_jumped_chr = true;

        uint num_variants = 0;
        uint num_variants_passed_non_de_novo = 0;
        uint num_variants_passed_de_novo = 0;

        mt19937 rng(12345678);

        while (true) {

            bool has_downstream_variant = in_vcf.getNextVariant(&cur_downstream_var);

            num_variants++;

            if (has_downstream_variant) {

                cur_downstream_var_chrom = cur_downstream_var->chrom();

                if (in_vcf.metaData().getContig(cur_downstream_var_chrom).type() != Contig::Type::Autosomal) {

                    continue;
                }

                if (cur_downstream_var_chrom == cur_eval_chr) {

                    downstream_var_jumped_chr = false;
                    downstream_var_start = cur_downstream_var->pos();

                } else {

                    downstream_var_jumped_chr = true;

                    if (!cur_eval_chr.empty()) {

                        downstream_var_start = in_vcf.metaData().getContig(cur_eval_chr).length() + 1;
                    }
                }

            } else {

                downstream_var_jumped_chr = true;
                downstream_var_start = in_vcf.metaData().getContig(cur_eval_chr).length() + 1;
            }

            if (cur_eval_var) {

                assert(downstream_var_start > 0);
                auto cur_eval_var_filters = cur_eval_var->filters();
                assert(cur_eval_var_filters.size() == 1);

                if (cur_eval_var_filters.front() == "PASS" and !(Auxiliaries::hasMissing(*cur_eval_var))) {

                    Allele ref_allele = cur_eval_var->ref();
                    int cur_eval_var_end = cur_eval_var->pos() + ref_allele.seq().size() - 1;

                    // Flank check
                    if (((static_cast<int>(cur_eval_var->pos()) - upstream_var_end - 1) >= min_flank_length) and ((downstream_var_start - cur_eval_var_end - 1) >= min_flank_length)) {

                        vector<uint> called_allele_idxs = Auxiliaries::getCalledAlleleIdxsSorted(*cur_eval_var, min_acp);

                        // Num called alleles and reference is called check
                        if (called_allele_idxs.size() == 2 and called_allele_idxs.front() == 0) {

                            uint val_alt_idx = called_allele_idxs.back();
                            Allele val_alt_allele = cur_eval_var->allele(val_alt_idx);

                            int max_allele_length = max(ref_allele.seq().size(), val_alt_allele.seq().size());

                            // No N's and fragment length check
                            if (!Auxiliaries::hasAmbiguous(*cur_eval_var) and ((max_allele_length + 2 * min_flank_length) <= max_fragment_length)) {

                                vector<ushort> het_all_idxs = {0, static_cast<ushort>(val_alt_idx)};

                                vector<string> hetzyg_max_post_sample_ids;
                                float hetzyg_max_post = 0;

                                for (auto sample_id : sample_ids) {

                                    auto cur_eval_var_sample = cur_eval_var->getSample(sample_id);

                                    if (cur_eval_var_sample.ploidy() == Sample::Ploidy::Diploid and cur_eval_var_sample.callStatus() == Sample::CallStatus::Complete) {

                                        auto gt_est = cur_eval_var_sample.genotypeEstimate();

                                        // If heterozygote ref/val_alt
                                        if ((find(gt_est.begin(), gt_est.end(), 0) != gt_est.end()) and (find(gt_est.begin(), gt_est.end(), val_alt_idx) != gt_est.end())) {

                                            auto gt_post = Auxiliaries::getGenotypePosterior(cur_eval_var_sample);
                                            assert(gt_post == Auxiliaries::getMaxGenotypePosterior(cur_eval_var_sample));
                                            
                                            assert(gt_post.second);

                                            if (floatCompare(gt_post.first, hetzyg_max_post)) {

                                                assert(hetzyg_max_post > 0);
                                                hetzyg_max_post_sample_ids.push_back(sample_id);

                                            } else if (gt_post.first > hetzyg_max_post) {

                                                hetzyg_max_post = gt_post.first;
                                                hetzyg_max_post_sample_ids.clear();
                                                hetzyg_max_post_sample_ids.push_back(sample_id);
                                            }
                                        }
                                    }
                                }

                                if (hetzyg_max_post >= min_gpp) {

                                    uniform_int_distribution<int> sample_idx_sampler(0, hetzyg_max_post_sample_ids.size() - 1);
                                    string hetzyg_max_post_sample_id = hetzyg_max_post_sample_ids.at(sample_idx_sampler(rng));

                                    // int max_flank = max_fragment_length - (max_allele_length + min_flank_length);
                                    int max_flank = (max_fragment_length - max_allele_length) / 2;
                                    assert(max_flank > 0);
                                    int upstream_flank_start = max(static_cast<int>(cur_eval_var->pos()) - max_flank, upstream_var_end + 1);
                                    assert(upstream_flank_start > 0);
                                    int downstream_flank_end = min(cur_eval_var_end + max_flank, downstream_var_start - 1);
                                    assert(downstream_flank_end > 0);

                                    auto genome_seq_find = genome_seqs.find(cur_eval_chr);
                                    assert(genome_seq_find != genome_seqs.end());

                                    string upstream_flank_seq = genome_seq_find->second->seq().substr(upstream_flank_start - 1, cur_eval_var->pos() - upstream_flank_start);
                                    assert(static_cast<int>(upstream_flank_seq.size()) <= max_flank);
                                    string downstream_flank_seq = genome_seq_find->second->seq().substr(cur_eval_var_end, downstream_flank_end - cur_eval_var_end);
                                    assert(static_cast<int>(downstream_flank_seq.size()) <= max_flank);

                                    string fastg_seq = upstream_flank_seq + "[" + ref_allele.seq() + "|" + val_alt_allele.seq() + "]" + downstream_flank_seq;

                                    auto total_allele_count = cur_eval_var->info().getValue<int>("AN");
                                    assert(total_allele_count.second);
                                    auto alt_allele_count = val_alt_allele.info().getValue<int>("AC");
                                    assert(alt_allele_count.second);

                                    auto allele_type = val_alt_allele.info().getValue<string>("AT");
                                    assert(allele_type.second);
                                    auto de_novo_info = cur_eval_var->info().getValue<string>("DNE");

                                    unordered_map<string,unique_ptr<ofstream> >::iterator out_file;

                                    if (de_novo_info.second) {

                                        num_variants_passed_de_novo++;
                                        out_file = output_files_denovo.find(allele_type.first);
                                        assert(out_file != output_files_denovo.end());

                                        // Assert that DNE info corresponds to var info
                                        auto dne_split = Utils::splitString(de_novo_info.first, ':');
                                        assert(dne_split.size() == 5);
                                        assert(dne_split.at(0) == hetzyg_max_post_sample_id);
                                        assert(dne_split.at(1) == "0");
                                        assert(dne_split.at(2) == to_string(val_alt_idx));
                                        assert(dne_split.at(3) == allele_type.first);

                                    } else {

                                        num_variants_passed_non_de_novo++;
                                        out_file = output_files.find(allele_type.first);
                                        assert(out_file != output_files.end());
                                    }

                                    Auxiliaries::rightTrimAllelePair(&ref_allele, &val_alt_allele);
                                    int max_trimmed_allele_length = max(ref_allele.seq().size(), val_alt_allele.seq().size());

                                    assert(max_allele_length > 0);
                                    assert(max_allele_length >= max_trimmed_allele_length);

                                    JoiningString val_var('\t');
                                    val_var.join(cur_eval_chr);
                                    val_var.join(to_string(cur_eval_var->pos()));
                                    val_var.join(to_string(max_allele_length));
                                    val_var.join(to_string(max_trimmed_allele_length));
                                    val_var.join(hetzyg_max_post_sample_id);
                                    val_var.join(to_string(hetzyg_max_post));
                                    val_var.join(to_string(upstream_flank_start));
                                    val_var.join(to_string(cur_eval_var->pos() - 1));
                                    val_var.join(to_string(cur_eval_var_end + 1));
                                    val_var.join(to_string(downstream_flank_end));
                                    val_var.join(fastg_seq);
                                    val_var.join(to_string(alt_allele_count.first));
                                    val_var.join(to_string(total_allele_count.first));

                                    *out_file->second << val_var.str() << endl;
                                }
                            }
                        }
                    }
                }
            }

            if (!has_downstream_variant) {

                break;
            }

            if (downstream_var_jumped_chr) {

                assert(!cur_downstream_var_chrom.empty());
                cur_eval_chr = cur_downstream_var_chrom;
                upstream_var_end = 0;

            } else {

                upstream_var_end = cur_eval_var->pos() + cur_eval_var->ref().seq().size() - 1;
            }

            if (cur_eval_var) delete cur_eval_var;
            cur_eval_var = cur_downstream_var;

            if ((num_variants % 100000) == 0) {

                std::cout << "[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants" << endl;
            }
        }

        for (auto & out_file : output_files) {

            out_file.second->close();
        }

        for (auto & out_file_denovo : output_files_denovo) {

            out_file_denovo.second->close();
        }

        for (auto & genome_seq : genome_seqs) {

            delete genome_seq.second;
        }

        cout << "[" << Utils::getLocalTime() << "] Completed BayesTyperUtils selectValVars\n" << endl;
        cout << "[" << Utils::getLocalTime() << "] " << num_variants << " variant(s) were parsed in total\n" << endl;
        cout << "[" << Utils::getLocalTime() << "] " << num_variants_passed_non_de_novo << " non de novo variant(s) were written to validation files\n" << endl;
        cout << "[" << Utils::getLocalTime() << "] " << num_variants_passed_de_novo << " de novo variant(s) were written to validation files\n" << endl;
    }
}
