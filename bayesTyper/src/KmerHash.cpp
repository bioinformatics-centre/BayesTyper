
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

using namespace std;

static const uint hybrid_hash_buckets_per_thread = 100000;

template<uchar kmer_size>
KmerHashHybrid<kmer_size>::KmerHashHybrid(const ushort num_threads, const ushort num_samples_in) : num_samples(num_samples_in) {

    min_kmer_count = 1;

	_hash = new ThreadedHybridHash<KmerCounts *, hash_root_size, kmer_size * 2 - hash_root_size>(hybrid_hash_buckets_per_thread * num_threads);
}


template<uchar kmer_size>
KmerHashHybrid<kmer_size>::~KmerHashHybrid() {

    cout << "\n[" << Utils::getLocalTime() << "] Deleting kmer hash ..." << endl;

    ulong kmer_count = 0;

    auto kit = _hash->begin();

    while (kit != _hash->end()) {

        if (*kit) {

            delete *kit;
            kmer_count++;
        }

        kit++;
    }

    delete _hash;

    cout << "[" << Utils::getLocalTime() << "] Deleted " << kmer_count << " kmers" << endl;
}


template<uchar kmer_size>
void KmerHashHybrid<kmer_size>::setMinKmerCount(const uchar min_kmer_count_in) {

	min_kmer_count = min_kmer_count_in;
}


template<uchar kmer_size>
ushort KmerHashHybrid<kmer_size>::getNumberOfSamples() {

    return num_samples;
}


template<uchar kmer_size>
vector<NegativeBinomialDistribution> KmerHashHybrid<kmer_size>::estimateGenomicCountDistributions(const vector<Sample> & samples) {

    cout << "\n[" << Utils::getLocalTime() << "] Tabulating intercluster kmer counts ..." << endl;

    vector<unordered_map<uint,ulong> > count_dists(samples.size());

    uint intercluster_counts = 0;
    auto kit = _hash->begin();

    while (kit != _hash->end()) {

        assert(*kit);

        if (!((*kit)->hasClusterOccurrence()) and !((*kit)->hasIncompleteClusterOccurrence())) {

            assert(!((*kit)->hasMulticlusterOccurrence()));
            assert(!((*kit)->hasMulticlusterGroupOccurrence()));

            if ((!(*kit)->hasMaxCount()) and (!(*kit)->hasDecoyOccurrence()) and ((*kit)->getInterclusterMultiplicity(Utils::Sex::Male) == 2) and ((*kit)->getInterclusterMultiplicity(Utils::Sex::Female) == 2)) {

                intercluster_counts++;

                for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

                    if ((*kit)->getCount(sample_idx) >= min_kmer_count) {

                        auto count_dists_ins_res = count_dists.at(sample_idx).emplace((*kit)->getCount(sample_idx), 0);
                        count_dists_ins_res.first->second++;
                    }
                }
            }

            delete *kit;
            *kit = nullptr;
    	}

        kit++;
    }

    cout << "[" << Utils::getLocalTime() << "] Finished tabulating kmer counts from " << intercluster_counts << " intercluster kmers:\n" << endl;

    assert(count_dists.size() == samples.size());

    vector<NegativeBinomialDistribution> genomic_count_distributions;
    genomic_count_distributions.reserve(count_dists.size());

    for (uint sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

        genomic_count_distributions.emplace_back(NegativeBinomialDistribution::methodOfMomentsEst(count_dists.at(sample_idx)));
        genomic_count_distributions.back().size(genomic_count_distributions.back().size()/2);

        cout << "[" << Utils::getLocalTime() << "] Estimated fixed negative binomial distribution for sample " << samples.at(sample_idx).name << " with mean " << genomic_count_distributions.back().mean() << " and variance " << genomic_count_distributions.back().var() << endl;
    }

    return genomic_count_distributions;
}

template<uchar kmer_size>
void KmerHashHybrid<kmer_size>::markKmers() {

    ulong total_count = 0;
    ulong unique_count = 0;
    ulong multicluster_count = 0;
    ulong intercluster_count = 0;
    ulong decoy_count = 0;
    ulong max_multiplicity_count = 0;
    ulong multicluster_group_count = 0;
    ulong constant_multiplicity_count = 0;
    ulong no_match_count = 0;

    auto kit = _hash->begin();

    while (kit != _hash->end()) {

        total_count++;

        if (*kit) {

            assert(!((*kit)->isEmpty()));

            if ((*kit)->hasClusterOccurrence()) {

                if ((*kit)->isExcluded()) {

                    if ((*kit)->hasDecoyOccurrence()) {

                        decoy_count++;

                    } else if ((*kit)->hasMaxMultiplicity()) {

                        max_multiplicity_count++;

                    } else if ((*kit)->hasMulticlusterGroupOccurrence()) {

                        multicluster_group_count++;
                    
                    } else {

                        assert(!((*kit)->hasMulticlusterOccurrence()));
                        assert((*kit)->hasConstantMultiplicity());

                        constant_multiplicity_count++;
                    } 

                } else if ((*kit)->hasMulticlusterOccurrence()) { 

                    assert(!((*kit)->hasMulticlusterGroupOccurrence()));

                    MultiClusterKmerCounts * multicluster_kmer_counts = new MultiClusterKmerCounts(num_samples, static_cast<SampleKmerCounts *>(*kit));
                    *kit = multicluster_kmer_counts;

                    multicluster_count++;

                } else if ((*kit)->getInterclusterMultiplicity(Utils::Sex::Male) > 0) {

                    intercluster_count++;

                } else {

                    unique_count++;
                }

            } else {

                assert((*kit)->hasIncompleteClusterOccurrence());
                no_match_count++;
            }

        } else {

            no_match_count++;
        }

        kit++;
    }

    assert(total_count == (unique_count + multicluster_count + intercluster_count + decoy_count + multicluster_group_count + max_multiplicity_count + constant_multiplicity_count + no_match_count));

    cout << "\n[" << Utils::getLocalTime() << "] Out of " << total_count << " unique kmers:\n" << endl;
    cout << "\t- " << unique_count << " have a unique match to a single variant cluster" << endl;
    cout << "\t- " << multicluster_count << " have a match to single variant cluster group and multiple internal variant clusters" << endl;
    cout << "\t- " << intercluster_count << " have a match to a single variant cluster and intercluster" << endl;
    
    cout << "\n\t- " << decoy_count << " have match to at least one variant cluster and has match to a decoy sequence (not used for inference)" << endl;
    cout << "\t- " << max_multiplicity_count << " have match to at least one variant cluster and has a maximum multiplicity higher than " << to_string(Utils::uchar_overflow) << " (not used for inference)" << endl;
    cout << "\t- " << multicluster_group_count << " have matches to multiple variant cluster groups (not used for inference)" << endl;
    cout << "\t- " << constant_multiplicity_count << " have match to a single variant cluster and has constant multiplicity across haplotype candidates (not used for inference)" << endl;
    
    cout << "\n\t- " << no_match_count << " have no match to a variant cluster\n" << endl;
}


template class KmerHashHybrid<31>;
template class KmerHashHybrid<39>;
template class KmerHashHybrid<47>;
template class KmerHashHybrid<55>;
template class KmerHashHybrid<63>;
