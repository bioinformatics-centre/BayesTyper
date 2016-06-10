
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
    };

    static pair<Stats::AlleleStats,bool> calcAlleleStats(Variant *, regex sample_id_incl_regex = regex(".+"), vector<AttributeFilter*> attribute_filters = vector<AttributeFilter*>());
    static pair<float,bool> calcInbreedingCoef(Variant *, regex sample_id_incl_regex = regex(".+"), vector<AttributeFilter*> attribute_filters = vector<AttributeFilter*>());
    static pair<vector<float>,float> calcAlleleCallProbAndQualFromAllelePosteriors(Variant *, regex sample_id_incl_regex = regex(".+"), vector<AttributeFilter*> attribute_filters = vector<AttributeFilter*>());


    class AlleleCounter {

        public:

            AlleleCounter(uint);
            virtual void addSample(Sample &, const vector<AttributeFilter*> &);
            uint alleleCountSum() const;
            const vector<uint> & alleleCounts() const;
            pair<vector<float>,bool> calcAlleleFreqs() const;

        protected:

            vector<uint> allele_counts;
            uint allele_count_sum;
    };

    class CompleteDiploidAlleleCounter : public AlleleCounter {

        public:

            CompleteDiploidAlleleCounter(uint);
            void addSample(Sample &, const vector<AttributeFilter*> &);

            pair<float,bool> calcNumExpectedHeterozygotes() const;
            uint numObservedHeterozygotes() const;

            uint num_heterozygotes;
    };

    template<typename Counter>
    static void applyCounterToSamples(Variant *, Counter *, const regex &, const vector<AttributeFilter*> &);
    static bool isAlleleFiltered(const AttributeSet &, const vector<AttributeFilter*> &);

};


#endif
