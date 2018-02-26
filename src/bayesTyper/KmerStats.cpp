
/*
KmerStats.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#include "KmerStats.hpp"
#include "KmerCounts.hpp"

using namespace std;

KmerStats::KmerStats() {

    count = 0;
    fraction = 0;
    mean = 0;
}

void KmerStats::reset() {

    count = 0;
    fraction = 0;
    mean = 0;
}

void KmerStats::addValue(const pair<float, bool> & value) {

    if (value.second) {

        count++;
        fraction += (static_cast<float>(!(Utils::floatCompare(value.first, 0))) - fraction) / count;
        mean += (value.first - mean) / count;
    }
}

uint KmerStats::getCount() const {

    return count;
}

pair<float, bool> KmerStats::getFraction() const {

    if (count == 0) {

        return make_pair(-1, false);

    } else {

        return make_pair(fraction, true);
    }
}

pair<float, bool> KmerStats::getMean() const {

    if (count == 0) {

        return make_pair(-1, false);

    } else {

        return make_pair(mean, true);
    }
}


AlleleKmerStats::AlleleKmerStats(const ushort num_alleles) {

    count_stats = vector<KmerStats>(num_alleles);
    fraction_stats = vector<KmerStats>(num_alleles);
    mean_stats = vector<KmerStats>(num_alleles);
}

void AlleleKmerStats::addKmerStats(const KmerStats & kmer_stats, const ushort allele_idx) {

    assert(allele_idx != Utils::ushort_overflow);

    count_stats.at(allele_idx).addValue(make_pair(kmer_stats.getCount(), true));
    fraction_stats.at(allele_idx).addValue(kmer_stats.getFraction());
    mean_stats.at(allele_idx).addValue(kmer_stats.getMean());

}



