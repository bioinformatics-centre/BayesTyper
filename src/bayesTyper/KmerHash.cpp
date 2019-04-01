
/*
KmerCountsHash.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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

#include "boost/iostreams/filtering_streambuf.hpp"
#include "boost/iostreams/filtering_stream.hpp"
#include "boost/iostreams/filter/gzip.hpp"

#include "KmerHash.hpp"
#include "Utils.hpp"
#include "KmerCounts.hpp"
#include "HybridHash.hpp"
#include "Sample.hpp"
#include "NegativeBinomialDistribution.hpp"
#include "Nucleotide.hpp"


using namespace std;

static const uint hash_root_size = pow(4,12);

static const uchar num_genomic_rate_gc_bias_bins = 1;
static const uchar max_multiplicity = Utils::uchar_overflow;


void writeSizeDistribution(const string & distribution_filename, const map<ulong, ulong> & size_distribution) {

    ofstream distribution_outfile(distribution_filename);

    if (!distribution_outfile.is_open()) {

        cerr << "\nERROR: Unable to write file " << distribution_filename << "\n" << endl;
        exit(1);
    }

    distribution_outfile << "Size\tCount" << endl;

    for (auto & size: size_distribution) {

        distribution_outfile << size.first << "\t" << size.second << endl;
    }

    distribution_outfile.close();
}

template<class T>
KmerHash<T>::KmerHash(const ulong expected_total_size, const ushort num_threads) {

    _hash = new ThreadedHybridHash<T, Utils::kmer_size * 2>(hash_root_size, expected_total_size, num_threads);
}

template<class T>
KmerHash<T>::~KmerHash() {

    delete _hash;
}

template<class T>
unique_lock<mutex> KmerHash<T>::getKmerLock(const bitset<Utils::kmer_size * 2> & kmer) {

    return move(_hash->lockKey(kmer));
}

template<class T>
pair<T *, bool> KmerHash<T>::addKmer(const bitset<Utils::kmer_size * 2> & kmer) {

    auto hash_insert = _hash->insert(kmer, T(0), true);
    assert(hash_insert.first != _hash->end());
   
    return make_pair(&((*hash_insert.first).second), hash_insert.second);
}

template<class T>
T * KmerHash<T>::findKmer(const bitset<Utils::kmer_size * 2> & kmer) {

    auto hash_it = _hash->find(kmer);

    if (hash_it != _hash->end()) {

        return &((*hash_it).second);

    } else {

        return nullptr;
    }
}

template<class T>
void KmerHash<T>::shuffle(const uint prng_seed) {

    _hash->shuffle(prng_seed);
}

template<class T>
ulong KmerHash<T>::size() {

    return _hash->size();
}

template<class T>
void KmerHash<T>::writeRootSizeDistribution(const string & output_prefix) {

    writeSizeDistribution(output_prefix + ".txt", _hash->getRootSizeDistribution());
}

template<class T>
ulong KmerHash<T>::writeKmersToFasta(const string & output_prefix, bool (*write_value)(T), const uint max_kmers) {

    ulong num_kmers_written = 0;

    ofstream kmers_outfile(output_prefix + ".fa.gz", ios::binary);

    if (!kmers_outfile.is_open()) {

        cerr << "\nERROR: Unable to write file " << output_prefix + ".fa.gz" << "\n" << endl;
        exit(1);
    }

    boost::iostreams::filtering_ostream kmers_outfile_fstream;

    kmers_outfile_fstream.push(boost::iostreams::gzip_compressor());
    kmers_outfile_fstream.push(boost::ref(kmers_outfile));

    assert(kmers_outfile_fstream.is_complete());    

    kmers_outfile_fstream << ">k" << Utils::kmer_size << endl;

    auto hash_it = _hash->begin();

    while (hash_it != _hash->end()) {

        if (write_value((*hash_it).second)) {

            kmers_outfile_fstream << Nucleotide::bitToNt<Utils::kmer_size>((*hash_it).first) << endl;
            num_kmers_written++;
        }

        if (num_kmers_written == max_kmers) {

            break;
        }

        hash_it++;
    }

    return num_kmers_written; 
}

template<class T>
ulong KmerHash<T>::addKmersToBloomFilter(KmerBloom<Utils::kmer_size> * kmer_bloom_filter, bool (*add_value)(T)) {

    ulong num_kmers_added = 0;

    auto hash_it = _hash->begin();

    while (hash_it != _hash->end()) {

        if (add_value((*hash_it).second)) {

            kmer_bloom_filter->addKmer((*hash_it).first);
            num_kmers_added++;
        }

        hash_it++;
    } 

    return num_kmers_added; 
}


template<uchar sample_bin>
ObservedKmerCountsHash<sample_bin>::ObservedKmerCountsHash(const ulong expected_total_size, const ushort num_threads) {

    _hash = new ThreadedHybridHash<ObservedKmerCounts<sample_bin>, Utils::kmer_size * 2>(hash_root_size, expected_total_size, num_threads);
}

template<uchar sample_bin>
ObservedKmerCountsHash<sample_bin>::~ObservedKmerCountsHash() {

    delete _hash;
}

template<uchar sample_bin>
unique_lock<mutex> ObservedKmerCountsHash<sample_bin>::getKmerLock(const bitset<Utils::kmer_size * 2> & kmer) {

    return move(_hash->lockKey(kmer));
}

template<uchar sample_bin>
pair<KmerCounts *, bool> ObservedKmerCountsHash<sample_bin>::addKmer(const bitset<Utils::kmer_size * 2> & kmer, const bool add_sorted) {

    auto hash_insert = _hash->insert(kmer, ObservedKmerCounts<sample_bin>(), add_sorted);
    assert(hash_insert.first != _hash->end());
   
    return make_pair(&((*hash_insert.first).second), hash_insert.second);
}

template<uchar sample_bin>
KmerCounts * ObservedKmerCountsHash<sample_bin>::findKmer(const bitset<Utils::kmer_size * 2> & kmer) {

    auto hash_it = _hash->find(kmer);

    if (hash_it != _hash->end()) {

        return &((*hash_it).second);

    } else {

        return nullptr;
    }
}

template<uchar sample_bin>
void ObservedKmerCountsHash<sample_bin>::sortKmers() {

    _hash->sort();
}

template<uchar sample_bin>
void ObservedKmerCountsHash<sample_bin>::writeRootSizeDistribution(const string & output_prefix) {

    writeSizeDistribution(output_prefix + ".txt", _hash->getRootSizeDistribution());
}

template<uchar sample_bin>
vector<vector<vector<KmerStats> > > ObservedKmerCountsHash<sample_bin>::calculateKmerStats(const vector<Sample> & samples) {

    vector<vector<vector<KmerStats> > > intercluster_kmer_stats(samples.size(), vector<vector<KmerStats> >(num_genomic_rate_gc_bias_bins, vector<KmerStats>(Utils::uchar_overflow + 1)));

    ulong total_count = 0;

    ulong unique_count = 0;
    ulong multicluster_count = 0;
    ulong decoy_count = 0;
    ulong max_multiplicity_count = 0;
    ulong multigroup_count = 0;
    ulong non_cluster_count = 0;
    
    auto hash_it = _hash->begin();

    while (hash_it != _hash->end()) {

        total_count++;

        if ((*hash_it).second.hasClusterOccurrence()) {

            assert(!(*hash_it).second.isParameter());

            if ((*hash_it).second.isExcluded()) {

                if ((*hash_it).second.hasDecoyOccurrence()) {

                    decoy_count++;

                } else if ((*hash_it).second.hasMaxMultiplicity()) {

                    max_multiplicity_count++;

                } else {

                    multigroup_count++;         
                    assert((*hash_it).second.hasMultigroupOccurrence());
                } 

            } else if ((*hash_it).second.hasMulticlusterOccurrence()) { 

                multicluster_count++;
                assert(!(*hash_it).second.hasMultigroupOccurrence());

            } else {

                unique_count++;
            }

        } else {

            non_cluster_count++;

            assert(!(*hash_it).second.hasMulticlusterOccurrence());
            assert(!(*hash_it).second.hasMultigroupOccurrence());

            if ((*hash_it).second.isParameter()) {

                const uchar bias_idx = Nucleotide::gcBiasBin<Utils::kmer_size>((*hash_it).first, num_genomic_rate_gc_bias_bins);

                for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {
            
                    intercluster_kmer_stats.at(sample_idx).at(bias_idx).at((*hash_it).second.getInterclusterMultiplicity(samples.at(sample_idx).gender)).addValue(make_pair((*hash_it).second.getSampleCount(sample_idx), true));
                }
            }
        }

        hash_it++;
    }

    assert(total_count == (unique_count + multicluster_count + decoy_count + multigroup_count + max_multiplicity_count + non_cluster_count));

    cout << "[" << Utils::getLocalTime() << "] Out of " << total_count << " kmers:\n" << endl;
    cout << "\t- " << unique_count << " have a match to a single variant cluster" << endl;
    cout << "\t- " << multicluster_count << " have a match to single variant cluster group and multiple variant clusters" << endl;
    
    cout << "\n\t- " << decoy_count << " have match to at least one variant cluster and has match to a decoy sequence (not used for inference)" << endl;
    cout << "\t- " << max_multiplicity_count << " have match to at least one variant cluster and has a maximum haploid multiplicity higher than " << to_string(Utils::bit7_overflow) << " (not used for inference)" << endl;
    cout << "\t- " << multigroup_count << " have matches to multiple variant cluster groups within or across inference units (not used for inference)" << endl;
    
    cout << "\n\t- " << non_cluster_count << " have no match to a variant cluster (includes parameter kmers)" << endl;
    
    return intercluster_kmer_stats;
}

template class KmerHash<bool>;
template class KmerHash<uchar>;

template class ObservedKmerCountsHash<3>;
template class ObservedKmerCountsHash<10>;
template class ObservedKmerCountsHash<20>;
template class ObservedKmerCountsHash<30>;

