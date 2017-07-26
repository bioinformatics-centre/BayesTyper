
/*
KmerStats.cpp - This file is part of BayesTyper (v0.9)


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

void KmerStats::addKmer(KmerCounts * const counts, const ushort sample_idx, const float multiplicity) {

    uchar kmer_count = 0;

    if (counts) {

        kmer_count = counts->getSampleCount(sample_idx);
    }

    assert(multiplicity > 0);

    count++;
    fraction += (static_cast<float>(kmer_count != 0) - fraction) / count;
    mean += ((kmer_count / multiplicity) - mean) / count;
}

float KmerStats::getCount() const {

    return count;
}

pair<float, bool> KmerStats::getFraction() const {

    if (Utils::floatCompare(count, 0)) {

        return make_pair(-1, false);

    } else {

        return make_pair(fraction, true);
    }
}

pair<float, bool> KmerStats::getMean() const {

    if (Utils::floatCompare(count, 0)) {

        return make_pair(-1, false);

    } else {

        return make_pair(mean, true);
    }
}


FixedKmerStats::FixedKmerStats() {

    count = -1;
    fraction = -1;
    mean = -1;
}

FixedKmerStats::FixedKmerStats(const float count_in, const float fraction_in, const float mean_in) {

    count = count_in;
    fraction = fraction_in;
    mean = mean_in;
}

float FixedKmerStats::getCount() const {

    return count;
}

float FixedKmerStats::getFraction() const {

    return fraction;
}

float FixedKmerStats::getMean() const {

    return mean;
}


MedianKmerStats::MedianKmerStats() {

    num_count_values = 0;
    num_fraction_values = 0;
    num_mean_values = 0;
}

void MedianKmerStats::addKmerStats(const KmerStats & kmer_stats) {

    auto count_values_it = count_values.insert(make_pair(kmer_stats.getCount(), 0));
    count_values_it.first->second++;
    num_count_values++;

    if (kmer_stats.getFraction().second) {

        assert(kmer_stats.getCount() > 0);

        auto fraction_values_it = fraction_values.insert(make_pair(kmer_stats.getFraction().first, 0));
        fraction_values_it.first->second++;
        num_fraction_values++;
    }

    if (kmer_stats.getMean().second) {

        assert(kmer_stats.getCount() > 0);

        auto mean_values_it = mean_values.insert(make_pair(kmer_stats.getMean().first, 0));
        mean_values_it.first->second++;
        num_mean_values++;
    }
}

pair<float, bool> MedianKmerStats::median(SortedLinearMap<float, uint, floatLess> & values, const uint num_values) {

    if (num_values == 0) {

        assert(values.empty());
        return make_pair(-1, false);

    } else {

        uint sum_num_values = 0;

        auto values_it = values.begin();
        assert(values_it != values.end());

        if ((num_values % 2) == 0) {

            const uint median_idx = static_cast<uint>(floor(static_cast<float>(num_values) / 2));

            while (values_it != values.end()) {

                sum_num_values += values_it->second;

                if (median_idx < sum_num_values) {

                    return make_pair(values_it->first, true);
                
                } else if (median_idx == sum_num_values) {

                    const float first_value = values_it->first;

                    values_it++;
                    assert(values_it != values.end());

                    return make_pair((first_value + values_it->first) / 2, true);
                }

                values_it++;
            }

        } else {

            const uint median_idx = static_cast<uint>(ceil(static_cast<float>(num_values) / 2));

            while (values_it != values.end()) {

                sum_num_values += values_it->second;

                if (median_idx <= sum_num_values) {

                    return make_pair(values_it->first, true);
                } 

                values_it++;
            }
        }

        assert(false);
        return make_pair(-1, false);
    }
}

pair<float, bool> MedianKmerStats::getMedianCount() {

    return median(count_values, num_count_values);
}

pair<float, bool> MedianKmerStats::getMedianFraction() {

    return median(fraction_values, num_fraction_values);
}

pair<float, bool> MedianKmerStats::getMedianMean() {
    
    return median(mean_values, num_mean_values);
}


VariantKmerStats::VariantKmerStats() {}

VariantKmerStats::VariantKmerStats(const ushort num_samples, const ushort num_alleles) {

    allele_kmer_stats = vector<vector<MedianKmerStats> >(num_samples, vector<MedianKmerStats>(num_alleles));
}

void VariantKmerStats::addAlleleKmerStats(const KmerStats & kmer_stats, const ushort allele_idx, const ushort sample_idx) {

    assert(allele_idx != Utils::ushort_overflow);
    allele_kmer_stats.at(sample_idx).at(allele_idx).addKmerStats(kmer_stats);
}



