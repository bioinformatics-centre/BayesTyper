
/*
KmerCounts.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#include <assert.h>
#include <vector>
#include <math.h>
#include <random>

#include "KmerCounts.hpp"
#include "Utils.hpp"

using namespace std;

KmerCounts::KmerCounts() {

    has_cluster_occ = false;
    has_multicluster_occ = false;
    has_multigroup_occ = false;

    has_decoy_occ = false;
    has_max_multiplicity = false;

    is_parameter = false;

    max_haploid_multiplicity = 0;

    female_intercluster_multiplicity = 0;
    male_intercluster_multiplicity = 0;
}


bool KmerCounts::hasClusterOccurrence() {

    return has_cluster_occ;
}

bool KmerCounts::hasMulticlusterOccurrence() {

    return has_multicluster_occ;
}

bool KmerCounts::hasMultigroupOccurrence() {

    return has_multigroup_occ;
}

bool KmerCounts::hasDecoyOccurrence() {

    return has_decoy_occ;
}

bool KmerCounts::hasMaxMultiplicity() {

    return has_max_multiplicity;
}

bool KmerCounts::isParameter() {

    return is_parameter;
}

void KmerCounts::isParameter(const bool is_parameter_in) {

    is_parameter = is_parameter_in;
}

bool KmerCounts::isExcluded() {

    return (has_decoy_occ or has_max_multiplicity or has_multigroup_occ);
}

void KmerCounts::addInterclusterMultiplicity(const bool is_decoy, const vector<Utils::Ploidy> & gender_ploidy) {

    max_haploid_multiplicity = updateMultiplicity(max_haploid_multiplicity, 1);

    if (max_haploid_multiplicity > Utils::bit7_overflow) {

        has_max_multiplicity = true;
    }

    if (is_decoy) {

        has_decoy_occ = true;

    } else {

        assert(gender_ploidy.size() == 2);

        female_intercluster_multiplicity = updateMultiplicity(female_intercluster_multiplicity, static_cast<uchar>(gender_ploidy.at(static_cast<uchar>(Utils::Gender::Female))));
        male_intercluster_multiplicity = updateMultiplicity(male_intercluster_multiplicity, static_cast<uchar>(gender_ploidy.at(static_cast<uchar>(Utils::Gender::Male))));
    }
}

uchar KmerCounts::getInterclusterMultiplicity(const Utils::Gender gender) {

    if (gender == Utils::Gender::Female) {

        return female_intercluster_multiplicity;
    
    } else {

        return male_intercluster_multiplicity;
    }
}

uchar KmerCounts::getMaxInterclusterMultiplicity() {

    return max(female_intercluster_multiplicity, male_intercluster_multiplicity);
}

void KmerCounts::addClusterMultiplicity(const uchar multiplicity, const bool is_multigroup_kmer) {

    if (has_cluster_occ) {

        has_multicluster_occ = true;
    }

    has_cluster_occ = true;

    if (is_multigroup_kmer) {

        has_multigroup_occ = true;
    }

    assert(has_multigroup_occ == is_multigroup_kmer);

    max_haploid_multiplicity = updateMultiplicity(max_haploid_multiplicity, multiplicity);     

    if (max_haploid_multiplicity > Utils::bit7_overflow) {

        has_max_multiplicity = true;
    }
}

uchar KmerCounts::updateMultiplicity(const uchar cur_multiplicity, const uchar in_multiplicity) {

    if ((Utils::uchar_overflow - cur_multiplicity) <= in_multiplicity) {

        return Utils::uchar_overflow;
    
    } else {

        return cur_multiplicity + in_multiplicity;
    } 
}


template <uchar sample_bin>
ObservedKmerCounts<sample_bin>::ObservedKmerCounts() : KmerCounts() {}

template <uchar sample_bin>
void ObservedKmerCounts<sample_bin>::addSampleCount(const ushort sample_idx, const uchar kmer_count) {

    if ((Utils::uchar_overflow - counts[sample_idx]) <= kmer_count) {

        counts[sample_idx] = Utils::uchar_overflow;    
    
    } else {

        counts[sample_idx] += kmer_count; 
    } 
}

template <uchar sample_bin>
uchar ObservedKmerCounts<sample_bin>::getSampleCount(const ushort sample_idx) {

    return counts[sample_idx];
}

template <uchar sample_bin>
void ObservedKmerCounts<sample_bin>::resetSampleMultiplicity() {

    for (ushort i = 0; i < sample_bin; i++) {

        multiplicities[i] = 0;
    }
}

template <uchar sample_bin>
uchar ObservedKmerCounts<sample_bin>::getSampleMultiplicity(const ushort sample_idx) {
    
    return multiplicities[sample_idx];
}

template <uchar sample_bin>
void ObservedKmerCounts<sample_bin>::reduceSampleMultiplicity(const ushort sample_idx, const uchar multiplicity) {
    
    assert(multiplicity <= multiplicities[sample_idx]);
    multiplicities[sample_idx] -= multiplicity;
}

template <uchar sample_bin>
void ObservedKmerCounts<sample_bin>::addSampleMultiplicity(const ushort sample_idx, const uchar multiplicity) {
    
    assert((multiplicity + multiplicities[sample_idx]) <= Utils::uchar_overflow);
    multiplicities[sample_idx] += multiplicity;
}

template class ObservedKmerCounts<3>;
template class ObservedKmerCounts<10>;
template class ObservedKmerCounts<20>;
template class ObservedKmerCounts<30>;



