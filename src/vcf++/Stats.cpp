
/*
Stats.cpp - This file is part of BayesTyper (v1.1)


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


bool Stats::isAlleleFiltered(Sample & cur_sample, const ushort allele_idx) {

    if (!(cur_sample.alleleInfo().empty())) {

        auto saf_value = cur_sample.alleleInfo().at(allele_idx).getValue<string>("SAF");

        if (saf_value.second) {

            assert((saf_value.first == "P") or (saf_value.first == "F"));

            if (saf_value.first == "F") {

                return true;
            }
        }
    }

    return false;
}

Stats::AlleleStats Stats::calcAlleleStats(Variant & variant, regex sample_id_incl_regex) {

    AlleleStats allele_stats(variant.numAlls());

    for (auto & sample_id: variant.sampleIds()) {

        if (regex_match(sample_id, sample_id_incl_regex)) {

            Sample & cur_sample = variant.getSample(sample_id);

            if (cur_sample.isInformative()) {

                auto cur_genotype = cur_sample.genotypeEstimate();

                for (auto & allele_idx: cur_genotype) {

                    if (!(isAlleleFiltered(cur_sample, allele_idx))) {

                        allele_stats.allele_count_sum++;
                        allele_stats.allele_counts.at(allele_idx)++;
                    }
                }
            }
        }
    }

    if (allele_stats.allele_count_sum > 0) {

        float allele_freq_cumsum = 0;

        for (ushort allele_idx = 0; allele_idx < variant.numAlls(); allele_idx++) {

            allele_stats.allele_freqs.at(allele_idx) = allele_stats.allele_counts.at(allele_idx)/static_cast<float>(allele_stats.allele_count_sum);
            allele_freq_cumsum += allele_stats.allele_freqs.at(allele_idx);
        }

        assert(floatCompare(allele_freq_cumsum, 1));
    }

    return allele_stats;
}

// From https://www.broadinstitute.org/gatk/guide/article?id=4732
Stats::InbreedingStats Stats::calcInbreedingStats(Variant & variant, regex sample_id_incl_regex) {

    vector<uint> allele_counts(variant.numAlls(), 0);

    uint num_eligible_samples = 0;
    uint num_obs_heterozygotes = 0;

    for (auto & sample_id : variant.sampleIds()) {

        Sample & cur_sample = variant.getSample(sample_id);

        if (regex_match(sample_id, sample_id_incl_regex) and cur_sample.ploidy() == Sample::Ploidy::Diploid and cur_sample.callStatus() == Sample::CallStatus::Complete) {

            auto sample_gt = cur_sample.genotypeEstimate();
            assert(sample_gt.size() == 2);

            if (!(isAlleleFiltered(cur_sample, sample_gt.front()) or isAlleleFiltered(cur_sample, sample_gt.back()))) {

                ++num_eligible_samples;

                allele_counts.at(sample_gt.front())++;
                allele_counts.at(sample_gt.back())++;

                if (sample_gt.front() != sample_gt.back()) {

                    ++num_obs_heterozygotes;
                }
            }
        }
    }

    InbreedingStats inbreeding_stats;
    inbreeding_stats.num_samples = num_eligible_samples;

    if (num_eligible_samples > 0) {

        float homozygote_prob = 0;

        for (auto & allele_count: allele_counts) {

            assert(allele_count <= (num_eligible_samples * 2));
            homozygote_prob += pow(allele_count/static_cast<float>(num_eligible_samples * 2), 2);
        }

        if (num_obs_heterozygotes != 0) {

            inbreeding_stats.is_fixed = false;

            float num_exp_heterozygotes = num_eligible_samples * (1 - homozygote_prob);
            assert(num_exp_heterozygotes > 0);

            inbreeding_stats.inbreeding_coef = 1 - (num_obs_heterozygotes/num_exp_heterozygotes);
        }
    }

    assert(Utils::floatCompare(inbreeding_stats.inbreeding_coef, -1) or (-1 < inbreeding_stats.inbreeding_coef));
    assert(Utils::floatCompare(inbreeding_stats.inbreeding_coef, 1) or (1 > inbreeding_stats.inbreeding_coef));

    return inbreeding_stats;
}


Stats::CallProbs Stats::calcCallProbs(Variant & variant, regex sample_id_incl_regex) {

    CallProbs call_probs(variant.numAlls());

    for (auto & sample_id: variant.sampleIds()) {

        if (regex_match(sample_id, sample_id_incl_regex)) {

            Sample & cur_sample = variant.getSample(sample_id);

            for (uint allele_idx = 0; allele_idx < cur_sample.alleleInfo().size(); allele_idx++) {

                if (!(isAlleleFiltered(cur_sample, allele_idx))) {

                    auto app_value = cur_sample.alleleInfo().at(allele_idx).getValue<float>("APP");

                    if (app_value.second) {

                        call_probs.allele_call_probs.at(allele_idx) = max(call_probs.allele_call_probs.at(allele_idx), app_value.first);
                    }
                }
            }
        }
    }

    float max_allele_call_prob = 0;

    for (uint allele_idx = 1; allele_idx < call_probs.allele_call_probs.size(); allele_idx++) {

        if (!(variant.allele(allele_idx).isMissing())) {

            max_allele_call_prob = max(max_allele_call_prob, call_probs.allele_call_probs.at(allele_idx));
        }
    }

    if (!(Utils::floatCompare(max_allele_call_prob, 0))) {

        if (Utils::floatCompare(max_allele_call_prob, 1)) {

            call_probs.variant_quality = 99;

        } else {

            assert(max_allele_call_prob > 0);
            assert(max_allele_call_prob < 1);

            call_probs.variant_quality = -10 * log10(1 - max_allele_call_prob);
        }
    }

    return call_probs;
}
