
/*
Stats.cpp - This file is part of BayesTyper (v0.9)


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


#include <numeric>
#include <algorithm>

#include "Stats.hpp"


Stats::AlleleCounter::AlleleCounter(uint num_alleles) {

    allele_counts = vector<uint>(num_alleles, 0);
    allele_count_sum = 0;
}

void Stats::AlleleCounter::addSample(Sample & sample, const vector<AttributeFilter*> & attribute_filters) {

    if (sample.isInformative()) {

        auto sample_gt = sample.genotypeEstimate();

        for (auto & gt_allele: sample_gt) {

            if (!(isAlleleFiltered(sample.alleleInfo().at(gt_allele), attribute_filters))) {

                allele_count_sum++;
                allele_counts.at(gt_allele)++;
            }
        }
    }
}

uint Stats::AlleleCounter::alleleCountSum() const {

    return allele_count_sum;
}

const vector<uint> & Stats::AlleleCounter::alleleCounts() const {

    return allele_counts;
}

pair<vector<float>,bool> Stats::AlleleCounter::calcAlleleFreqs() const {

    if (allele_count_sum > 0) {

        vector<float> allele_freqs;
        allele_freqs.reserve(allele_counts.size());
        float allele_freq_cum = 0;

        for (uint allele_count : allele_counts) {

            float allele_freq = float(allele_count)/float(allele_count_sum);
            allele_freqs.push_back(allele_freq);
            allele_freq_cum += allele_freq;
        }

        assert(floatCompare(allele_freq_cum, 1));

        return make_pair(allele_freqs, true);

    } else {

        return make_pair(vector<float>(allele_counts.size(), 0), false);
    }
}

Stats::CompleteDiploidAlleleCounter::CompleteDiploidAlleleCounter(uint num_alleles) : AlleleCounter(num_alleles) {

    num_heterozygotes = 0;    
}

void Stats::CompleteDiploidAlleleCounter::addSample(Sample & sample, const vector<AttributeFilter*> & attribute_filters) {

    if (sample.ploidy() == Sample::Ploidy::Diploid and sample.callStatus() == Sample::CallStatus::Complete) {

        auto sample_gt = sample.genotypeEstimate();
        assert(sample_gt.size() == 2);

        if (!(isAlleleFiltered(sample.alleleInfo().at(sample_gt.front()), attribute_filters) or isAlleleFiltered(sample.alleleInfo().at(sample_gt.back()), attribute_filters))) {

            allele_count_sum += 2;

            allele_counts.at(sample_gt.front())++;
            allele_counts.at(sample_gt.back())++;

            if (sample_gt.front() != sample_gt.back()) {

                num_heterozygotes++;
            }
        }
    }
}

pair<float,bool> Stats::CompleteDiploidAlleleCounter::calcNumExpectedHeterozygotes() const {

    if (allele_count_sum > 0) {

        uint num_alleles = allele_counts.size();
        auto allele_freqs = calcAlleleFreqs();
        assert(allele_freqs.second);
        assert(allele_freqs.first.size() == num_alleles);

        float exp_hetero_frac = 0;

        for (uint i = 0; i < num_alleles; i++) {

            for (uint j = i + 1; j < num_alleles; j++) {

                exp_hetero_frac += 2 * allele_freqs.first.at(i) * allele_freqs.first.at(j);
            }
        }

        float num_samples = float(allele_count_sum)/2; // Due to diploid
        float exp_num_hetero = exp_hetero_frac * num_samples;

        if (Utils::floatCompare(exp_num_hetero, 0)) {

            return make_pair(0, false);

        } else {

            return make_pair(exp_num_hetero, true);
        }

    } else {

        return make_pair(0, false);
    }
}

uint Stats::CompleteDiploidAlleleCounter::numObservedHeterozygotes() const {

    return num_heterozygotes;
}

template<typename Counter>
void Stats::applyCounterToSamples(Variant * variant, Counter * counter, const regex & sample_id_incl_regex, const vector<AttributeFilter*> & attribute_filters) {

    for (auto & sample_id: variant->sampleIds()) {

        if (regex_match(sample_id, sample_id_incl_regex)) {

            counter->addSample(variant->getSample(sample_id), attribute_filters);
        }
    }
}

bool Stats::isAlleleFiltered(const AttributeSet & attribute_set, const vector<AttributeFilter*> & attribute_filters) {

    auto aff_value = attribute_set.getValue<string>("AFF");

    if (aff_value.second) {

        assert((aff_value.first == "P") or (aff_value.first == "F"));

        if (aff_value.first == "F") {

            return true;
        }
    }

    for (auto attribute_filter: attribute_filters) {

        if (!(attribute_filter->pass(attribute_set))) {

            return true;
        }
    }

    return false; 
}

pair<Stats::AlleleStats,bool> Stats::calcAlleleStats(Variant * variant, regex sample_id_incl_regex, vector<AttributeFilter*> attribute_filters) {

    AlleleCounter allele_counter(variant->numAlls());
    applyCounterToSamples(variant, &allele_counter, sample_id_incl_regex, attribute_filters);

    AlleleStats allele_stats;

    allele_stats.allele_count_sum = allele_counter.alleleCountSum();
    allele_stats.allele_counts = allele_counter.alleleCounts();
    
    auto allele_freqs = allele_counter.calcAlleleFreqs();
    allele_stats.allele_freqs = allele_freqs.first;

    assert(allele_stats.allele_counts.size() == variant->numAlls());
    assert(allele_stats.allele_freqs.size() == variant->numAlls());

    if (allele_counter.alleleCountSum() > 0) {

        assert(allele_freqs.second);

        return make_pair(allele_stats, true);

    } else {

        assert(!allele_freqs.second);

        return make_pair(allele_stats, false);
    }
}

// From https://www.broadinstitute.org/gatk/guide/article?id=4732, generalised to > biallelic
pair<float,bool> Stats::calcInbreedingCoef(Variant * variant, regex sample_id_incl_regex, vector<AttributeFilter*> attribute_filters) {

    CompleteDiploidAlleleCounter allele_counter(variant->numAlls());
    applyCounterToSamples(variant, &allele_counter, sample_id_incl_regex, attribute_filters);

    auto num_exp_heterozygotes = allele_counter.calcNumExpectedHeterozygotes();
    uint num_obs_heterozygotes = allele_counter.numObservedHeterozygotes();

    if (num_exp_heterozygotes.second) {

        return make_pair(1-(float(num_obs_heterozygotes)/num_exp_heterozygotes.first), true);

    } else {

        return make_pair(0, false);
    }
}

pair<vector<float>,float> Stats::calcAlleleCallProbAndQualFromAllelePosteriors(Variant * cur_var, regex sample_id_incl_regex, vector<AttributeFilter*> attribute_filters) {

    vector<float> mapp_call_prob(cur_var->numAlls(), 0);
    float max_call_prob = 0;

    auto sample_ids = cur_var->sampleIds();

    for (uint all_idx = 0; all_idx < cur_var->numAlls(); all_idx++) {

        double not_called_prob = 0;
        bool has_perfect = false;

        for (auto & id : sample_ids) {

            if (regex_match(id, sample_id_incl_regex)) {
        
                Sample & cur_sample = cur_var->getSample(id); 

                if (!(cur_sample.isInformative())) {

                    continue;
                } 

                if (!(isAlleleFiltered(cur_sample.alleleInfo().at(all_idx), attribute_filters))) {

                    auto mapp_value = cur_sample.alleleInfo().at(all_idx).getValue<float>("MAP");
                    assert(mapp_value.second);

                    if (!(Utils::floatCompare(mapp_value.first, 1))) {

                        not_called_prob += log(1 - mapp_value.first);

                    } else {

                        has_perfect = true;
                        break;
                    }
                }
            }
        }

        if (has_perfect) {

            mapp_call_prob.at(all_idx) = 1;

        } else {

            mapp_call_prob.at(all_idx) = 1 - exp(not_called_prob);
        }

        if (all_idx > 0 and !cur_var->allele(all_idx).isMissing()) {

            max_call_prob = max(max_call_prob, mapp_call_prob.at(all_idx));
        }
    }

    float qual;

    if (Utils::floatCompare(max_call_prob, 1)) {

        qual = 999;

    } else if (Utils::floatCompare(max_call_prob, 0)) {

        qual = 0;

    } else {

        assert(max_call_prob > 0);
        assert(max_call_prob < 1);

        qual = -10*log10(1 - max_call_prob);
    }

    return make_pair(mapp_call_prob, qual);
}

