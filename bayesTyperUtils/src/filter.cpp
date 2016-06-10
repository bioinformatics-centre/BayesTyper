
/*
filter.cpp - This file is part of BayesTyper (v0.9)


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


#include <algorithm>

#include "vcf++/VcfFile.hpp"
#include "vcf++/Stats.hpp"
#include "vcf++/Auxiliaries.hpp"

#include "filter.hpp"

namespace Filter {

    void filter(const string & in_vcf_filename, const string & out_vcf_prefix, const vector<pair<AttributeFilter*, AttributeFilter*> > & attribute_filters, const float min_gpp, const bool keep_filtered_variants, const bool filter_dependencies) {

        cout << "\n[" << Utils::getLocalTime() << "] Running BayesTyperUtils (" << BTU_VERSION << ") filter with " << attribute_filters.size() + 1 << " filter(s) ...\n" << endl;

        GenotypedVcfFileReader in_vcf(in_vcf_filename, true);
        VcfMetaData output_meta_data = in_vcf.metaData();

        assert(output_meta_data.filterDescriptors().emplace("ISD", Attribute::Descriptor("ISD", "Variant filtered due to insufficient data")).second);

        VcfFileWriter out_vcf(out_vcf_prefix + ".vcf", output_meta_data, true);

        Variant * cur_var;
        auto sample_ids = output_meta_data.sampleIds();

        string cur_chromosome = "";

        unordered_map<string, uint> sample_filtered_alleles;

        for (auto & sample_id: sample_ids) {

            assert(sample_filtered_alleles.emplace(sample_id, 0).second);
        }

        ulong num_variants = 0;
        ulong num_variants_passed = 0;
        ulong num_variants_filtered_uv = 0;
        ulong num_variants_filtered_isd = 0;

        ulong num_genotypes = 0;
        ulong num_genotypes_passed = 0;
        ulong num_genotypes_filtered_pre = 0;
        ulong num_genotypes_filtered_par = 0;
        ulong num_genotypes_filtered_ful = 0;

        while (in_vcf.getNextVariant(&cur_var)) {

            num_variants++;

            assert(cur_var->filters().size() == 1);

            if (cur_var->filters().front() != "PASS") {

                assert(cur_var->filters().front() == "UV");
                num_variants_filtered_uv++;

                num_genotypes += cur_var->numSamples();
                num_genotypes_filtered_pre += cur_var->numSamples();

                if (keep_filtered_variants) {

                    out_vcf.write(cur_var);
                }

                delete cur_var;
                continue;
            }

            if (cur_var->chrom() != cur_chromosome) {

                cur_chromosome = cur_var->chrom();

                for (auto & sample: sample_filtered_alleles) {

                    sample.second = 0;
                }
            }

            unordered_set<uint> excl_alleles;
            auto cur_var_ea = cur_var->info().getValue<string>("EA");

            if (cur_var_ea.second) {

                auto cur_var_ea_split = Utils::splitString(cur_var_ea.first, ',');

                for (auto & excl_allele_idx: cur_var_ea_split) {

                    assert(excl_alleles.insert(stoi(excl_allele_idx)).second);
                }
            }

            for (auto & sample_id: sample_ids) {

                num_genotypes++;

                Sample * cur_sample = &(cur_var->getSample(sample_id));
                auto cur_genotype_estimate = cur_sample->genotypeEstimate();

                bool is_pre_filtered = false;

                if (cur_sample->callStatus() == Sample::CallStatus::Missing) {

                    is_pre_filtered = true;

                    assert(cur_genotype_estimate.empty());
                    num_genotypes_filtered_pre++;
                }

                for (auto & attribute_set: cur_sample->alleleInfo()) {

                    bool is_filtered = false;

                    for (auto attribute_filter: attribute_filters) {

                        if (!(attribute_filter.first->pass(attribute_set) or attribute_filter.second->pass(attribute_set))) {

                            is_filtered = true;
                            break;
                        }
                    }

                    if (is_filtered) {

                        assert(!attribute_set.setValue<string>("AFF", "F"));
                    }             
                }

                vector<ushort> new_genotype_estimate;
                new_genotype_estimate.reserve(cur_genotype_estimate.size());

                for (auto allele_idx: cur_genotype_estimate) {

                    auto aff_value = cur_sample->alleleInfo().at(allele_idx).getValue<string>("AFF");
                    assert(aff_value.second);

                    if (aff_value.first == "P") {

                        new_genotype_estimate.push_back(allele_idx);
                    }
                }

                assert(new_genotype_estimate.size() <= cur_genotype_estimate.size());

                if (!(new_genotype_estimate.empty())) {

                    auto max_gpp = Auxiliaries::getMaxGenotypePosterior(*cur_sample);
                    assert(max_gpp.second);

                    if (cur_sample->callStatus() == Sample::CallStatus::Complete) {

                        assert(max_gpp == Auxiliaries::getGenotypePosterior(*cur_sample));
                    }

                    if (max_gpp.first < min_gpp) {

                        new_genotype_estimate.clear();
                    } 
                }

                auto sample_filtered_alleles_it = sample_filtered_alleles.find(sample_id);
                assert(sample_filtered_alleles_it != sample_filtered_alleles.end());

                if (cur_var->pos() <= sample_filtered_alleles_it->second) {

                    assert(Auxiliaries::hasMissing(*cur_var));
                    assert(filter_dependencies);

                    new_genotype_estimate.clear();
                }

                if (new_genotype_estimate.size() < cur_genotype_estimate.size()) {

                    cur_sample->newGenotypeEstimate(new_genotype_estimate);
                }

                if (filter_dependencies) {

                    for (uint alt_idx = 0; alt_idx < cur_var->numAlts(); alt_idx++) {

                        auto aff_value = cur_sample->alleleInfo().at(alt_idx + 1).getValue<string>("AFF");
                        assert(aff_value.second);

                        if (aff_value.first == "P") {

                            continue;
                        }

                        if (cur_var->ref().seq().size() == 1) {

                            continue;
                        }

                        if (cur_var->alt(alt_idx).isMissing()) {

                            continue;
                        }

                        if (excl_alleles.find(alt_idx + 1) != excl_alleles.end()) {

                            continue;
                        }

                        Allele cur_ref_allele = cur_var->ref();
                        Allele cur_alt_allele = cur_var->alt(alt_idx);

                        Auxiliaries::rightTrimAllelePair(&cur_ref_allele, &cur_alt_allele);

                        assert(!(cur_ref_allele.seq().empty()));
                        assert(!(cur_alt_allele.seq().empty()));

                        sample_filtered_alleles_it->second = max(sample_filtered_alleles_it->second, static_cast<uint>(cur_var->pos() + cur_ref_allele.seq().size() - 1));
                    }
                }

                if (!is_pre_filtered) {

                    if (cur_sample->callStatus() == Sample::CallStatus::Complete) {

                        num_genotypes_passed++;
                    
                    } else if (cur_sample->callStatus() == Sample::CallStatus::Partial) {

                        num_genotypes_filtered_par++;

                    } else {

                        assert(cur_sample->callStatus() == Sample::CallStatus::Missing);
                        num_genotypes_filtered_ful++;
                    }
                }
            }

            auto allele_call_prob_post = Stats::calcAlleleCallProbAndQualFromAllelePosteriors(cur_var);
            cur_var->setQual(make_pair(allele_call_prob_post.second, true));

            auto allele_stats = Stats::calcAlleleStats(cur_var);

            assert(!cur_var->info().setValue<int>("AN", allele_stats.first.allele_count_sum));

            for (uint allele_idx = 0; allele_idx < cur_var->numAlls(); allele_idx++) {

                assert(!cur_var->allele(allele_idx).info().setValue<float>("ACP", allele_call_prob_post.first.at(allele_idx)));

                if (allele_idx > 0) {

                    assert(!cur_var->allele(allele_idx).info().setValue<int>("AC", allele_stats.first.allele_counts.at(allele_idx)));
                    assert(!cur_var->allele(allele_idx).info().setValue<float>("AF", allele_stats.first.allele_freqs.at(allele_idx)));
                }
            }

            if (allele_stats.first.allele_count_sum > 0) {

                num_variants_passed++;
                out_vcf.write(cur_var);

            } else {

                num_variants_filtered_isd++;
                cur_var->setFilters(vector<string>(1,"ISD"));

                if (keep_filtered_variants) {

                    out_vcf.write(cur_var);
                }
            }

            delete cur_var;

            if ((num_variants % 100000) == 0) {

                std::cout << "[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants" << endl;
            }
        }

        cout << "\n[" << Utils::getLocalTime() << "] Completed filtering of " << num_variants << " variants";
        cout <<  " with writing of filtered variants";

        if (keep_filtered_variants) {

            cout << " enabled" << endl;

        } else {

            cout << " disabled" << endl;
        }

        cout << "\n[" << Utils::getLocalTime() << "] " << num_variants_passed << " variants passed filters" << endl;
        cout << "[" << Utils::getLocalTime() << "] " << num_variants - num_variants_passed << " variants failed:\n" << endl;
        cout << "\t- " << num_variants_filtered_uv << " variants filtered due to lack of support for variant type (UV)" << endl;
        cout << "\t- " << num_variants_filtered_isd << " variants filtered due to insufficient data (ISD)" << endl;

        assert(num_variants == (num_variants_passed + num_variants_filtered_uv + num_variants_filtered_isd));

        cout << "\n[" << Utils::getLocalTime() << "] " << num_genotypes_passed << " genotypes passed filters" << endl;
        cout << "[" << Utils::getLocalTime() << "] " << num_genotypes - num_genotypes_passed << " genotypes failed:\n" << endl;
        cout << "\t- " << num_genotypes_filtered_pre << " genotypes were pre-filtered (fully filterd in input set)" << endl;
        cout << "\t- " << num_genotypes_filtered_par << " genotypes were partial filtered" << endl;
        cout << "\t- " << num_genotypes_filtered_ful << " genotypes were fully filtered" << endl;

        cout << endl;

        assert(num_genotypes == (num_genotypes_passed + num_genotypes_filtered_pre + num_genotypes_filtered_par + num_genotypes_filtered_ful));
    }
}
