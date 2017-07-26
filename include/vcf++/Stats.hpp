
/*
Stats.hpp - This file is part of BayesTyper (v0.9)


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


#ifndef STATS
#define STATS

#include <regex>

#include "AttributeFilter.hpp"
#include "SampleAlleleAttributeFilter.hpp"
#include "Variant.hpp"
#include "Sample.hpp"
#include "Utils.hpp"

class Stats {

public:

    struct AlleleStats {

        uint allele_count_sum;
        vector<uint> allele_counts;
        vector<float> allele_freqs;

        AlleleStats(const ushort num_alleles) {

            allele_count_sum = 0;
            allele_counts = vector<uint>(num_alleles, 0);
            allele_freqs = vector<float>(num_alleles, 0);
        }
    };

    struct InbreedingStats {

        uint num_samples;
        float inbreeding_coef;
        bool is_fixed;

        InbreedingStats() {

            num_samples = 0;
            inbreeding_coef = 0;
            is_fixed = true;
        }
    };


    struct CallProbs {

        float variant_quality;
        vector<float> allele_call_probs;

        CallProbs(const ushort num_alleles) {

            variant_quality = 0;
            allele_call_probs = vector<float>(num_alleles, 0);
        }
    };

    static Stats::AlleleStats calcAlleleStats(Variant &, regex sample_id_incl_regex = regex(".+"));
    static Stats::InbreedingStats calcInbreedingStats(Variant &, regex sample_id_incl_regex = regex(".+"));
    static Stats::CallProbs calcCallProbs(Variant &, regex sample_id_incl_regex = regex(".+"));

    static bool isAlleleFiltered(Sample &, const ushort);

};


#endif
