
/*
KmerHash.cpp - This file is part of BayesTyper (v0.9)


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
#include <map>
#include <string>
#include <math.h>

#include "boost/math/special_functions/gamma.hpp"

#include "KmerHash.hpp"
#include "Utils.hpp"
#include "KmerCounts.hpp"
#include "HybridHash.hpp"
#include "Sample.hpp"
#include "NegativeBinomialDistribution.hpp"
#include "Sequence.hpp"

using namespace std;


template<uchar kmer_size, uchar sample_bin>
HybridKmerHash<kmer_size, sample_bin>::HybridKmerHash(const ushort num_samples_in, const ushort num_threads) : num_samples(num_samples_in) {

    _hash = new ThreadedHybridHash<ObservedKmerCounts<sample_bin>, hash_root_size, kmer_size * 2 - hash_root_size>(num_threads);
}

template<uchar kmer_size, uchar sample_bin>
HybridKmerHash<kmer_size, sample_bin>::~HybridKmerHash() {

    delete _hash;
}

template<uchar kmer_size, uchar sample_bin>
unique_lock<mutex> HybridKmerHash<kmer_size, sample_bin>::getKmerLock(const bitset<kmer_size * 2> & kmer) {

    return move(_hash->lockKey(kmer));
}

template<uchar kmer_size, uchar sample_bin>
pair<KmerCounts *, bool> HybridKmerHash<kmer_size, sample_bin>::addKmer(const bitset<kmer_size * 2> & kmer, const bool add_sorted) {

    auto hash_insert = _hash->insert(kmer, ObservedKmerCounts<sample_bin>(), add_sorted);
    assert(hash_insert.first != _hash->end());
   
    return make_pair(&((*hash_insert.first).second), hash_insert.second);
}

template<uchar kmer_size, uchar sample_bin>
KmerCounts * HybridKmerHash<kmer_size, sample_bin>::findKmer(const bitset<kmer_size * 2> & kmer) {

    auto hash_it = _hash->find(kmer);

    if (hash_it != _hash->end()) {

        return &((*hash_it).second);

    } else {

        return nullptr;
    }
}

template<uchar kmer_size, uchar sample_bin>
void HybridKmerHash<kmer_size, sample_bin>::sortKmers() {

    _hash->sort();
}

template<uchar kmer_size, uchar sample_bin>
vector<vector<vector<ulong> > > HybridKmerHash<kmer_size, sample_bin>::calculateKmerStats(const uchar num_genomic_rate_gc_bias_bins) {

    vector<vector<vector<ulong> > > diploid_kmer_counts(num_samples, vector<vector<ulong> >(num_genomic_rate_gc_bias_bins, vector<ulong>(Utils::uchar_overflow + 1, 0)));

    ulong total_count = 0;
    ulong unique_count = 0;
    ulong multicluster_count = 0;
    ulong decoy_count = 0;
    ulong max_multiplicity_count = 0;
    ulong multicluster_group_count = 0;
    ulong constant_multiplicity_count = 0;
    ulong no_match_count = 0;

    auto kit = _hash->begin();

    while (kit != _hash->end()) {

        total_count++;

        if ((*kit).second.hasClusterOccurrence()) {

            if ((*kit).second.isExcluded()) {

                if ((*kit).second.hasDecoyOccurrence()) {

                    decoy_count++;

                } else if ((*kit).second.hasMaxMultiplicity()) {

                    max_multiplicity_count++;

                } else if ((*kit).second.hasMulticlusterGroupOccurrence()) {

                    multicluster_group_count++;
                
                } else {

                    assert((*kit).second.hasConstantMultiplicity());
                    constant_multiplicity_count++;
                } 

            } else if ((*kit).second.hasMulticlusterOccurrence()) { 

                assert(!((*kit).second.hasMulticlusterGroupOccurrence()));
                multicluster_count++;

            } else {

                unique_count++;
            }

        } else {

            no_match_count++;

            assert(!((*kit).second.hasMulticlusterOccurrence()));
            assert(!((*kit).second.hasMulticlusterGroupOccurrence()));

            if (!((*kit).second.hasDecoyOccurrence()) and ((*kit).second.getInterclusterMultiplicity(Utils::Gender::Male) == 2) and ((*kit).second.getInterclusterMultiplicity(Utils::Gender::Female) == 2)) {

                for (ushort sample_idx = 0; sample_idx < num_samples; sample_idx++) {

                    diploid_kmer_counts.at(sample_idx).at(Sequence::gcBiasBin<kmer_size>((*kit).first, num_genomic_rate_gc_bias_bins)).at((*kit).second.getSampleCount(sample_idx))++;
                }
            }
        }

        kit++;
    }

    assert(total_count == (unique_count + multicluster_count + decoy_count + multicluster_group_count + max_multiplicity_count + constant_multiplicity_count + no_match_count));

    cout << "[" << Utils::getLocalTime() << "] Out of " << total_count << " unique kmers:\n" << endl;
    cout << "\t- " << unique_count << " have a unique match to a single variant cluster" << endl;
    cout << "\t- " << multicluster_count << " have a match to single variant cluster group and multiple internal variant clusters" << endl;
    
    cout << "\n\t- " << decoy_count << " have match to at least one variant cluster and has match to a decoy sequence (not used for inference)" << endl;
    cout << "\t- " << max_multiplicity_count << " have match to at least one variant cluster and has a maximum haploid multiplicity higher than " << to_string(Utils::bit7_overflow) << " (not used for inference)" << endl;
    cout << "\t- " << multicluster_group_count << " have matches to multiple variant cluster groups (not used for inference)" << endl;
    cout << "\t- " << constant_multiplicity_count << " have match to at least one variant cluster and has constant multiplicity across all haplotype candidates (not used for inference)" << endl;
    
    cout << "\n\t- " << no_match_count << " have no match to a variant cluster\n" << endl;

    return diploid_kmer_counts;
}


template class BasicKmerHash<31>;
template class BasicKmerHash<39>;
template class BasicKmerHash<47>;
template class BasicKmerHash<55>;
template class BasicKmerHash<63>;


template class HybridKmerHash<31,1>;
template class HybridKmerHash<39,1>;
template class HybridKmerHash<47,1>;
template class HybridKmerHash<55,1>;
template class HybridKmerHash<63,1>;

template class HybridKmerHash<31,3>;
template class HybridKmerHash<39,3>;
template class HybridKmerHash<47,3>;
template class HybridKmerHash<55,3>;
template class HybridKmerHash<63,3>;

template class HybridKmerHash<31,10>;
template class HybridKmerHash<39,10>;
template class HybridKmerHash<47,10>;
template class HybridKmerHash<55,10>;
template class HybridKmerHash<63,10>;

template class HybridKmerHash<31,20>;
template class HybridKmerHash<39,20>;
template class HybridKmerHash<47,20>;
template class HybridKmerHash<55,20>;
template class HybridKmerHash<63,20>;

template class HybridKmerHash<31,30>;
template class HybridKmerHash<39,30>;
template class HybridKmerHash<47,30>;
template class HybridKmerHash<55,30>;
template class HybridKmerHash<63,30>;
