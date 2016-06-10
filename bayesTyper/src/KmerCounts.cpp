
/*
KmerCounts.cpp - This file is part of BayesTyper (v0.9)


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

    has_max_count = false;
    has_cluster_occ = false;
    has_multicluster_occ = false;
    has_incomplete_cluster_occ = false;

    has_decoy_occ = false;
    has_max_multiplicity = false;
    has_multicluster_group_occ = false;
    has_constant_multiplicity = false;
}


bool KmerCounts::hasMaxCount() {

    return has_max_count;
}

bool KmerCounts::hasClusterOccurrence() {

    return has_cluster_occ;
}

bool KmerCounts::hasMulticlusterOccurrence() {

    return has_multicluster_occ;
}

bool KmerCounts::hasIncompleteClusterOccurrence() {

    return has_incomplete_cluster_occ;
}

bool KmerCounts::hasDecoyOccurrence() {

    return has_decoy_occ;
}

bool KmerCounts::hasMaxMultiplicity() {

    return has_max_multiplicity;
}

bool KmerCounts::hasMulticlusterGroupOccurrence() {

    return has_multicluster_group_occ;
}

bool KmerCounts::hasConstantMultiplicity() {

    return has_constant_multiplicity;
}

bool KmerCounts::isExcluded() {

    return (has_decoy_occ or has_max_multiplicity or has_multicluster_group_occ or (!has_multicluster_occ and has_constant_multiplicity));
}


SampleKmerCounts::SampleKmerCounts(const ushort num_samples) : KmerCounts() {

    counts = new uchar[num_samples];

    for (auto i = 0; i < num_samples; i++) {
        
        counts[i] = 0;
    }

    max_cluster_multiplicity = 0; 
    male_intercluster_multiplicity = 0;
    female_intercluster_multiplicity = 0;

    variant_cluster_group_index = Utils::uint_overflow;
}

SampleKmerCounts::~SampleKmerCounts() {

    delete[] counts;
}

void SampleKmerCounts::addCount(const ushort sample_idx, const uchar kmer_count) {

    if ((Utils::uchar_overflow - counts[sample_idx]) <= kmer_count) {

        has_max_count = true;
        counts[sample_idx] = Utils::uchar_overflow;    
    
    } else {

        counts[sample_idx] += kmer_count; 
    } 
}

uchar SampleKmerCounts::addMultiplicity(const uchar cur_multiplicity, const uchar in_multiplicity) {

    if ((Utils::uchar_overflow - cur_multiplicity) <= in_multiplicity) {

        has_max_multiplicity = true;
        return Utils::uchar_overflow;
    
    } else {

        return cur_multiplicity + in_multiplicity;
    } 
}

void SampleKmerCounts::addClusterMultiplicity(const uchar multiplicity, const bool has_constant_multiplicity_in, const uint variant_cluster_group_index_in) {

    if (multiplicity == 0) {

        has_incomplete_cluster_occ = true;
    
    } else {

        if (has_cluster_occ) {

            has_multicluster_occ = true;
        }

        has_cluster_occ = true;

        if (has_constant_multiplicity_in) {

            has_constant_multiplicity = has_constant_multiplicity_in;
        }

        if (multiplicity > Utils::bit7_overflow) {

            max_cluster_multiplicity = addMultiplicity(max_cluster_multiplicity, Utils::uchar_overflow);     

        } else {

            max_cluster_multiplicity = addMultiplicity(max_cluster_multiplicity, 2 * multiplicity);     
        } 

        if (!has_max_multiplicity) {

            if ((max_cluster_multiplicity + max(male_intercluster_multiplicity, female_intercluster_multiplicity)) > Utils::uchar_overflow) {

                has_max_multiplicity = true;
            }
        }

        if (variant_cluster_group_index == Utils::uint_overflow) {

            variant_cluster_group_index = variant_cluster_group_index_in;

        } else {

            assert(has_multicluster_occ);

            if (variant_cluster_group_index != variant_cluster_group_index_in) {

                has_multicluster_group_occ = true;
            }         
        }
    } 
}

void SampleKmerCounts::addInterclusterMultiplicity(const Utils::ChromosomeClass chromosome_class) {

    if (chromosome_class == Utils::ChromosomeClass::Autosomal) {

        male_intercluster_multiplicity = addMultiplicity(male_intercluster_multiplicity, 2);
        female_intercluster_multiplicity = addMultiplicity(female_intercluster_multiplicity, 2);
    
    } else if (chromosome_class == Utils::ChromosomeClass::X) {

        male_intercluster_multiplicity = addMultiplicity(male_intercluster_multiplicity, 1);
        female_intercluster_multiplicity = addMultiplicity(female_intercluster_multiplicity, 2);

    } else if (chromosome_class == Utils::ChromosomeClass::Y) {

        male_intercluster_multiplicity = addMultiplicity(male_intercluster_multiplicity, 1);
    
    } else {

        assert(chromosome_class == Utils::ChromosomeClass::Decoy);  
        has_decoy_occ = true;
    }
}

uchar SampleKmerCounts::getCount(const ushort sample_idx) {

    return counts[sample_idx];
}

uchar SampleKmerCounts::getInterclusterMultiplicity(const Utils::Sex sex) {

    if (sex == Utils::Sex::Male) {

        return male_intercluster_multiplicity;
    
    } else {

        return female_intercluster_multiplicity;
    }
}

bool SampleKmerCounts::isMulti() {

    return false;
}

bool SampleKmerCounts::isEmpty() {

    return false;
}


MultiClusterKmerCounts::MultiClusterKmerCounts(const ushort num_samples, SampleKmerCounts * sample_kmer_counts_in) : KmerCounts(), sample_kmer_counts(sample_kmer_counts_in) {

    assert(sample_kmer_counts);

    multiplicities = new uchar[num_samples];
    indices = new uint[num_samples];

    for (auto i = 0; i < num_samples; i++) {
        
        multiplicities[i] = 0;
        indices[i] = Utils::ushort_overflow;
    }

    has_max_count = sample_kmer_counts->hasMaxCount();
    has_cluster_occ = sample_kmer_counts->hasClusterOccurrence();
    has_multicluster_occ = sample_kmer_counts->hasMulticlusterOccurrence();
    has_incomplete_cluster_occ = sample_kmer_counts->hasIncompleteClusterOccurrence();

    has_decoy_occ = sample_kmer_counts->hasDecoyOccurrence();
    has_max_multiplicity = sample_kmer_counts->hasMaxMultiplicity();
    has_multicluster_group_occ = sample_kmer_counts->hasMulticlusterGroupOccurrence();
    has_constant_multiplicity = sample_kmer_counts->hasConstantMultiplicity();

    assert(has_cluster_occ);
    assert(has_multicluster_occ);

    assert(!has_decoy_occ);
    assert(!has_max_multiplicity);
    assert(!has_multicluster_group_occ);
}

MultiClusterKmerCounts::~MultiClusterKmerCounts() {

    delete sample_kmer_counts;
    delete[] multiplicities;
    delete[] indices;
}

void MultiClusterKmerCounts::reset(const ushort num_samples) {
    
    for (auto i = 0; i < num_samples; i++) {
        
        multiplicities[i] = 0;
    }
}

uchar MultiClusterKmerCounts::getMulticlusterMultiplicity(const ushort sample_idx) {

    if (sample_kmer_counts->getCount(sample_idx) == 0) {

        return 0;
    
    } else {

        return multiplicities[sample_idx];
    }
}

void MultiClusterKmerCounts::reduceMulticlusterMultiplicity(const ushort sample_idx, const uchar multiplicity) {
    
    assert(multiplicity <= multiplicities[sample_idx]);
    multiplicities[sample_idx] -= multiplicity;
}

void MultiClusterKmerCounts::addMulticlusterMultiplicity(const ushort sample_idx, const uchar multiplicity) {
    
    assert((multiplicity + multiplicities[sample_idx]) <= Utils::uchar_overflow);
    multiplicities[sample_idx] += multiplicity;
}

uchar MultiClusterKmerCounts::getCount(const ushort sample_idx) {

    return sample_kmer_counts->getCount(sample_idx);
}

uchar MultiClusterKmerCounts::getInterclusterMultiplicity(const Utils::Sex sex) {
    
    return sample_kmer_counts->getInterclusterMultiplicity(sex);
}

void MultiClusterKmerCounts::setIndex(const ushort sample_idx, const uint index_in) {

    indices[sample_idx] = index_in;
}

uint MultiClusterKmerCounts::getIndex(const ushort sample_idx) {

    return indices[sample_idx];
}

bool MultiClusterKmerCounts::isMulti() {

    return true;
}

bool MultiClusterKmerCounts::isEmpty() {

    return false;
}


EmptyKmerCounts::EmptyKmerCounts() : KmerCounts() {}

uchar EmptyKmerCounts::getCount(const ushort) {

    return 0;
}

uchar EmptyKmerCounts::getInterclusterMultiplicity(const Utils::Sex sex) {

    return 0;
}

bool EmptyKmerCounts::isMulti() {

    return false;
}

bool EmptyKmerCounts::isEmpty() {

    return true;
}



