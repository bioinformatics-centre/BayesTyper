
/*
VariantClusterGraph.cpp - This file is part of BayesTyper (v0.9)


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


#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <bitset>
#include <list>
#include <mutex>
#include <math.h> 
#include <random>
#include <algorithm>

#include "boost/graph/adjacency_list.hpp"
#include "../Eigen/Dense"

#include "VariantClusterGraph.hpp"
#include "Utils.hpp"
#include "VariantCluster.hpp"
#include "KmerHash.hpp"
#include "PerfectSet.hpp"
#include "VariantClusterGraphKmerPath.hpp"
#include "Kmer.hpp"
#include "KmerCounts.hpp"
#include "VariantClusterHaplotypes.hpp"
#include "VariantKmerStats.hpp"
#include "VariantInfo.hpp"
#include "VariantClusterGraphVertex.hpp"
#include "KmerInfo.hpp"

static const ulong max_small_complexity = pow(4, Utils::small_kmer_size) * 0.001;
static const ushort max_recursion_depth = 100;

static const ushort haplotype_candidate_nbest_multiplicity = 3;
static const ushort num_dummy_haplotype_candidates = 24;


VariantClusterGraph::VariantClusterGraph(VariantCluster * variant_cluster, string & chromosome_sequence, const uchar kmer_size) {

    has_ambiguous_nucleotide = false;
    has_complex_region = false;
    has_redundant_sequence = false;
    has_non_unique_kmer = false;
    has_excluded_kmer = false;

    assert(variant_cluster->variants.size() < Utils::ushort_overflow);
    
    assert(variant_cluster_info.empty());
    variant_cluster_info.reserve(variant_cluster->variants.size());

    map<uint, pair<vector<vertex_t>, vector<ushort> > > added_vertices;
    unordered_set<ushort> reference_variant_ids;

    auto lit = variant_cluster->variants.begin();
    auto cit = chromosome_sequence.begin();

    vertex_t cur_vertex = boost::add_vertex(graph);    
    vertex_t prev_vertex = cur_vertex;

    addVertex(cur_vertex, vector<StringItPair>(1, StringItPair(cit + lit->first - (kmer_size - 1), cit + lit->first)), vector<uint>(), make_pair(Utils::ushort_overflow, Utils::ushort_overflow), reference_variant_ids, kmer_size);

    assert(added_vertices.insert({lit->first, make_pair(vector<vertex_t>(1, cur_vertex), vector<ushort>())}).second);
    pair<uint, vertex_t> prev_new_position(lit->first, cur_vertex);

    uint cur_last_position = 0;
    uint next_position = 0;

    ushort variant_counter = 0;

    while (lit != variant_cluster->variants.end()) {

        // assert(variant_ids.insert(lit->first).second);

        variant_cluster_info.emplace_back(VariantInfo());
        variant_cluster_info.back().variant_id = lit->second.id;
        variant_cluster_info.back().variant_type = lit->second.type;
        variant_cluster_info.back().num_alleles = lit->second.alternative_alleles.size() + 1 + lit->second.excluded_alternative_alleles.size();
        variant_cluster_info.back().excluded_alternative_alleles = lit->second.excluded_alternative_alleles;
        variant_cluster_info.back().has_dependency = lit->second.has_dependency;

        for (auto & ait: lit->second.alternative_alleles) {

            assert(count(lit->second.excluded_alternative_alleles.begin(), lit->second.excluded_alternative_alleles.end(), ait.first) == 0);

            vertex_t next_vertex = boost::add_vertex(graph);
            addVertex(next_vertex, vector<StringItPair>(1, StringItPair(ait.second.second.begin(), ait.second.second.end())), vector<uint>(), make_pair(variant_counter, ait.first), reference_variant_ids, kmer_size);

            if (prev_new_position.first == lit->first) { 

                assert(boost::add_edge(prev_new_position.second, next_vertex, graph).second);

            } else {

                assert(boost::add_edge(cur_vertex, next_vertex, graph).second);
            }

            auto added_vertices_insert = added_vertices.insert({lit->first + ait.second.first, make_pair(vector<vertex_t>(), vector<ushort>())});
            added_vertices_insert.first->second.first.push_back(next_vertex);
        }

        if (prev_new_position.first != lit->first) {

            prev_new_position.first = lit->first;
            prev_new_position.second = cur_vertex;
        }

        auto avit = added_vertices.find(lit->first + lit->second.max_reference_length);
        assert(avit != added_vertices.end());
        avit->second.second.push_back(variant_counter);

        reference_variant_ids.insert(variant_counter);

        lit++;

        auto added_vertices_it = added_vertices.begin();
        bool more_edges = true;
        bool last_variant = false;

        if (lit != variant_cluster->variants.end()) {

            next_position = lit->first;

        } else {

            next_position = Utils::uint_overflow;
            last_variant = true;
        }

        while (more_edges) {

            auto cur_position = added_vertices_it->first;
            auto next_vertices = added_vertices_it->second.first;

            for (auto &variant_id: added_vertices_it->second.second) {

                reference_variant_ids.erase(variant_id);
            }

            added_vertices.erase(added_vertices_it);

            if (added_vertices.empty()) {

                more_edges = false;

                if (last_variant) {

                    cur_last_position = cur_position + kmer_size - 1;

                } else {

                    cur_last_position = next_position;
                }                   

            } else {

                added_vertices_it = added_vertices.begin();
                cur_last_position = added_vertices_it->first;

                if ((!last_variant) and (cur_last_position > next_position)) {

                    more_edges = false;
                    cur_last_position = next_position;
                }
            }
            
            vector<StringItPair> contained_vertices;
            contained_vertices.reserve(variant_cluster->contained_clusters.size());
            
            vector<uint> nested_variant_cluster_indices;

            auto contained_cluster_it = variant_cluster->contained_clusters.begin();
            uint prev_contained_edge = Utils::uint_overflow;

            while ((contained_cluster_it != variant_cluster->contained_clusters.end()) and (contained_cluster_it->left_flank < cur_last_position)) {

                assert(cur_position <= contained_cluster_it->left_flank);                
                assert(contained_cluster_it->right_flank <= (cur_last_position - kmer_size));

                if (prev_contained_edge < Utils::uint_overflow) {

                    nested_variant_cluster_indices.emplace_back(prev_contained_edge);
                }

                contained_vertices.emplace_back(cit + cur_position, cit + contained_cluster_it->left_flank);

                prev_contained_edge = contained_cluster_it->cluster_idx;
                cur_position = contained_cluster_it->right_flank + 1;
                variant_cluster->contained_clusters.erase(contained_cluster_it++);
            }
 
            if (prev_contained_edge < Utils::uint_overflow) {

                nested_variant_cluster_indices.emplace_back(prev_contained_edge);
            }

            contained_vertices.emplace_back(cit + cur_position, cit + cur_last_position);

            cur_vertex = boost::add_vertex(graph);
            bool is_reference_allele = false;

            for (auto &vit: next_vertices) {

                if (vit == prev_vertex) {

                    is_reference_allele = true;
                    assert(next_vertices.size() == 1);
                }

                assert(boost::add_edge(vit, cur_vertex, graph).second);
            }

            if (is_reference_allele) {

                addVertex(cur_vertex, contained_vertices, nested_variant_cluster_indices, pair<ushort, ushort>(variant_counter, 0), reference_variant_ids, kmer_size);

            } else {

                addVertex(cur_vertex, contained_vertices, nested_variant_cluster_indices, pair<ushort, ushort>(Utils::ushort_overflow, Utils::ushort_overflow), reference_variant_ids, kmer_size);
            }

            auto added_vertices_insert = added_vertices.insert({cur_last_position, make_pair(vector<vertex_t>(), vector<ushort>())});
            added_vertices_insert.first->second.first.push_back(cur_vertex);
        }

        variant_counter++;
        prev_vertex = cur_vertex; 
    }

    assert(variant_cluster->variants.size() == variant_counter);

    assert(added_vertices.size() == 1);
    assert(added_vertices.begin()->second.first.size() == 1);
    assert(added_vertices.begin()->second.second.empty());
}


VariantClusterGraph::~VariantClusterGraph() {

    auto vit = boost::vertices(graph);
    assert(vit.first != vit.second);

    while (vit.first != vit.second) {

        delete graph[*vit.first];
        vit.first++;
    }
}


void VariantClusterGraph::addVertex(vertex_t cur_vertex, vector<StringItPair> vertex_sequnces, vector<uint> nested_variant_cluster_indices, pair<ushort, ushort> alt_allele_id, unordered_set<ushort> & ref_allele_ids, const uchar kmer_size) {

    vector<vector<bool> > cur_sequences;

    assert(!(vertex_sequnces.empty()));

    if (vertex_sequnces.size() == 1) {

        addSequence(&cur_sequences, vertex_sequnces.front().first, vertex_sequnces.front().second);
    
    } else {

        addSequence(&cur_sequences, vertex_sequnces.front().first, vertex_sequnces.front().second);
        
        for (uint i = 1; i < (vertex_sequnces.size() - 1); i++) {
    
            addSequence(&cur_sequences, vertex_sequnces.at(i).first, vertex_sequnces.at(i).second);
        }

        addSequence(&cur_sequences, vertex_sequnces.back().first, vertex_sequnces.back().second);
    }

    assert(!(cur_sequences.empty()));

    if (nested_variant_cluster_indices.empty()) {

        graph[cur_vertex] = new VariantClusterGraphVertex();

    } else {

        assert(nested_variant_cluster_indices.size() <= cur_sequences.size());

        graph[cur_vertex] = new ComplexVariantClusterGraphVertex();
        static_cast<ComplexVariantClusterGraphVertex * >(graph[cur_vertex])->nested_variant_cluster_indices = nested_variant_cluster_indices;
    }

    graph[cur_vertex]->sequences = move(cur_sequences);

    if (alt_allele_id.first != Utils::ushort_overflow) {

        assert(alt_allele_id.second != Utils::ushort_overflow);

        graph[cur_vertex]->variant_ids.emplace_back(VariantClusterGraphVertex::VariantId(alt_allele_id, true));
    }

    for (auto &ref_id: ref_allele_ids) {

        if (ref_id != alt_allele_id.first) {

            graph[cur_vertex]->variant_ids.emplace_back(VariantClusterGraphVertex::VariantId(pair<ushort, ushort>(ref_id, 0), false));
        }
    }
}


void VariantClusterGraph::addSequence(vector<vector<bool> > * cur_sequences, string::const_iterator sit, string::const_iterator eit) {

    cur_sequences->emplace_back(vector<bool>());
    cur_sequences->back().reserve((eit - sit) * 2);

    bitset<2> nt_bitmer;
    bool last_not_N = true;

    while (sit != eit) {

        if (!(Utils::DNAtoBit(*sit, &nt_bitmer))) {

            has_ambiguous_nucleotide = true;

            if (last_not_N) {

                cur_sequences->back().shrink_to_fit();
                assert(cur_sequences->back().size()%2 == 0);

                cur_sequences->emplace_back(vector<bool>());
                cur_sequences->back().reserve((eit - sit) * 2);
            }

            last_not_N = false;
            sit++;;

            continue;
        }
        
        cur_sequences->back().push_back(nt_bitmer[0]);
        cur_sequences->back().push_back(nt_bitmer[1]);

        last_not_N = true;
        sit++;
    }

    assert(cur_sequences->back().size()%2 == 0);
}


ulong VariantClusterGraph::countSmallmers(Utils::SmallmerSet * smallmer_set) {

    auto vit = boost::vertices(graph);

    auto in_eit = boost::in_edges(*vit.first, graph);
    auto out_eit = boost::out_edges(*vit.first, graph);

    assert(*vit.first == 0);
    assert(in_eit.first == in_eit.second);
    assert(out_eit.first != out_eit.second);

    ushort max_alternative_alleles = Utils::small_kmer_size + 1;
    ulong cur_max_small_complexity = 0;

    do {

        max_alternative_alleles--;
        vit = boost::vertices(graph);

        cur_max_small_complexity = 0;

        while (vit.first != vit.second) {

            cur_max_small_complexity = max(cur_max_small_complexity, calculateVertexSmallmerComplexity(*vit.first, 0, 0, max_alternative_alleles));

            if (cur_max_small_complexity > max_small_complexity) {

                break;
            }

            vit.first++;
        }
    
    } while (cur_max_small_complexity > max_small_complexity);

    if (max_alternative_alleles == Utils::small_kmer_size) {

        max_alternative_alleles = Utils::ushort_overflow;
    
    } else {

        has_complex_region = true;
    }

    vit = boost::vertices(graph);
    set<vertex_t> not_visited_vertices;

    while (vit.first != vit.second) {

        assert(not_visited_vertices.insert(*vit.first).second);
        vit.first++;
    }

    assert(*not_visited_vertices.begin() == 0);

    ulong num_unique_smallmers = 0;

    while (!(not_visited_vertices.empty())) {

        KmerPair<Utils::small_kmer_size> cur_pair;
        num_unique_smallmers += countVertexSmallmers(smallmer_set, cur_pair, Utils::uint_overflow, *not_visited_vertices.begin(), &not_visited_vertices, 0, max_alternative_alleles, max_recursion_depth);
    }

    return num_unique_smallmers;
}


ulong VariantClusterGraph::calculateVertexSmallmerComplexity(const vertex_t cur_vertex, uint cur_small_mer_size, uchar num_alternative_alleles, const ushort max_alternative_alleles) {

    if (!(graph[cur_vertex]->variant_ids.empty())) {

        if ((graph[cur_vertex]->variant_ids.front().is_allele_vertex) and (graph[cur_vertex]->variant_ids.front().allele_id.second > 0)) {

            num_alternative_alleles++;
        }
    }

    if (num_alternative_alleles == max_alternative_alleles) {

        return 0;
    }

    if (cur_small_mer_size == 0) {

        cur_small_mer_size = 1;

    } else {

        for (uint seq_idx = 0; seq_idx < graph[cur_vertex]->sequences.size(); seq_idx++) {

            cur_small_mer_size += graph[cur_vertex]->sequences.at(seq_idx).size()/2;
        }
    }

    if (cur_small_mer_size >= Utils::small_kmer_size) {

        return 1;
    }

    ulong num_small_mers = 0;

    auto out_eit = boost::out_edges(cur_vertex, graph);

    while (out_eit.first != out_eit.second) {

        auto target_vertex = boost::target(*out_eit.first, graph);
        assert(target_vertex > cur_vertex);

        num_small_mers += calculateVertexSmallmerComplexity(target_vertex, cur_small_mer_size, num_alternative_alleles, max_alternative_alleles);

        out_eit.first++;
    }

    return num_small_mers;
}


ulong VariantClusterGraph::countVertexSmallmers(Utils::SmallmerSet * smallmer_set, KmerPair<Utils::small_kmer_size> cur_pair, uint remaining_prefix, const vertex_t cur_vertex, set<vertex_t> * unvisited_vertices, ushort num_alternative_alleles, const ushort max_alternative_alleles, ushort cur_depth) {

    if (!(graph[cur_vertex]->variant_ids.empty())) {

        if ((graph[cur_vertex]->variant_ids.front().is_allele_vertex) and (graph[cur_vertex]->variant_ids.front().allele_id.second > 0)) {

            num_alternative_alleles++;
        }
    }

    if (num_alternative_alleles == max_alternative_alleles) {

        return 0;
    }

    if (remaining_prefix >= Utils::small_kmer_size) {

        unvisited_vertices->erase(cur_vertex);
    }

    ulong num_unique_smallmers = 0;

    bitset<2> nt_bits;

    for (uint seq_idx = 0; seq_idx < graph[cur_vertex]->sequences.size(); seq_idx++) {

        if (seq_idx > 0) {

            cur_pair.reset();
        }

        auto sit = graph[cur_vertex]->sequences.at(seq_idx).cbegin();

        while (sit != graph[cur_vertex]->sequences.at(seq_idx).cend()) {

            nt_bits.set(0, *sit);
            sit++;

            nt_bits.set(1, *sit);
            sit++;

           if (cur_pair.shift(nt_bits)) { 

                num_unique_smallmers += smallmer_set->insert(cur_pair.getForwardKmer());
                num_unique_smallmers += smallmer_set->insert(cur_pair.getReverseComplementKmer());
            } 

            remaining_prefix--;

            if (remaining_prefix == 0) {

                return num_unique_smallmers;
            } 
        }
    }

    auto out_eit = boost::out_edges(cur_vertex, graph);

    while (out_eit.first != out_eit.second) {
                    
        if (remaining_prefix < Utils::small_kmer_size) {

            num_unique_smallmers += countVertexSmallmers(smallmer_set, cur_pair, remaining_prefix, boost::target(*out_eit.first, graph), unvisited_vertices, num_alternative_alleles, max_alternative_alleles, cur_depth);

        } else if (unvisited_vertices->find(boost::target(*out_eit.first, graph)) == unvisited_vertices->end()) {

            num_unique_smallmers += countVertexSmallmers(smallmer_set, cur_pair, Utils::small_kmer_size - 1, boost::target(*out_eit.first, graph), unvisited_vertices, num_alternative_alleles, max_alternative_alleles, cur_depth);

        } else if (cur_depth == 0) {

            num_unique_smallmers += countVertexSmallmers(smallmer_set, cur_pair, Utils::small_kmer_size - 1, boost::target(*out_eit.first, graph), unvisited_vertices, num_alternative_alleles, max_alternative_alleles, cur_depth);

        } else {

            num_unique_smallmers += countVertexSmallmers(smallmer_set, cur_pair, Utils::uint_overflow, boost::target(*out_eit.first, graph), unvisited_vertices, num_alternative_alleles, max_alternative_alleles, cur_depth - 1);
        }

        out_eit.first++;
    }

    return num_unique_smallmers;
}


const vector<VariantInfo> & VariantClusterGraph::getInfo() {

    return variant_cluster_info;  
}


bool VariantClusterGraph::hasAmbiguousNucleotide() {

    return has_ambiguous_nucleotide;
}


bool VariantClusterGraph::hasComplexRegion() {

    return has_complex_region;
}


bool VariantClusterGraph::hasRedundantSequence() {

    return has_redundant_sequence;
}


bool VariantClusterGraph::hasNonUniqueKmer() {

    return has_non_unique_kmer;
}


bool VariantClusterGraph::hasExcludedKmer() {

    return has_excluded_kmer;
}


// void VariantClusterGraph::printInfoAndGraph() {

//     cout << "\n" << endl;

//     auto vit = vertices(graph);

//     while (vit.first != vit.second) {

//         cout << *vit.first << " : [";

//         for (auto &id: graph[*vit.first]->variant_ids) {

//             cout << "(" << id.allele_id.first << ", " << id.allele_id.second << ", " << id.is_allele_vertex << "), ";
//         }

//         cout << "] [";

//         for (auto &sequence: graph[*vit.first]->sequences) {

//             cout << sequence << "(" << sequence.size() << "), ";
//         }

//         cout << "] [";

//         if (graph[*vit.first]->isComplex()) {

//             for (auto &id: static_cast<ComplexVariantClusterGraphVertex * >(graph[*vit.first])->nested_variant_cluster_indices) {

//                 cout << id << ", ";
//             }
//         }

//         cout << "]" << endl;

//         vit.first++;
//     }

//     cout << "\n" << endl;

//     auto eit = edges(graph);

//     while (eit.first != eit.second) {

//         cout << boost::source(*eit.first, graph) << "->" << boost::target(*eit.first, graph) << endl;
//         eit.first++;
//     }

//     cout << endl;
// }


template <uchar kmer_size>
VariantClusterGraphLockedKmerSize<kmer_size>::VariantClusterGraphLockedKmerSize(VariantCluster * variant_cluster, string & chromosome_sequence) : VariantClusterGraph(variant_cluster, chromosome_sequence, kmer_size) {
}


template <uchar kmer_size>
ulong VariantClusterGraphLockedKmerSize<kmer_size>::countKmers(KmerHash * kmer_hash, const uint prng_seed, const ushort num_haplotype_candidates_per_sample, const uint variant_cluster_group_idx, double * variant_cluster_group_complexity) {

    KmerHashHybrid<kmer_size> * kmer_hash_hybrid = static_cast<KmerHashHybrid<kmer_size> * >(kmer_hash);
    
    mt19937 prng = mt19937(prng_seed);
    
    vector<KmerMultiplicities *> path_kmer_multiplicities;
    path_kmer_multiplicities.reserve(variant_cluster_info.size() * 2 * num_haplotype_candidates_per_sample * (kmer_hash_hybrid->getNumberOfSamples() + 1));

    vector<vector<VariantClusterGraphKmerPath<kmer_size> *> > * best_paths_per_sample; 

    best_paths_per_sample = findBestPathsPerSample<VariantClusterGraphKmerPath<kmer_size> >(kmer_hash_hybrid, kmer_hash_hybrid->getNumberOfSamples(), &path_kmer_multiplicities, &prng, num_haplotype_candidates_per_sample);   

    assert(best_paths_per_sample->size() == static_cast<uint>(kmer_hash_hybrid->getNumberOfSamples() + 1));

    vector<VariantClusterGraphKmerPath<kmer_size> *> collapsed_best_paths = collapsePaths(best_paths_per_sample, num_haplotype_candidates_per_sample * (kmer_hash_hybrid->getNumberOfSamples() + 1));
    delete best_paths_per_sample;

    assert(!(collapsed_best_paths.empty()));

    unordered_set<KmerMultiplicities *> visited_kmer_vertices;
    visited_kmer_vertices.reserve(path_kmer_multiplicities.size());

    KmerMultiplicitiesIndex kmer_multiplicities_index = createKmerMultiplicitiesIndex(&collapsed_best_paths, &visited_kmer_vertices, true);

    double number_of_paths = collapsed_best_paths.size();
    *variant_cluster_group_complexity += kmer_multiplicities_index.size() * (((number_of_paths - 1) * number_of_paths) / 2 + number_of_paths);

    ulong num_counted_sample_kmers = 0;

    for (auto &kit: kmer_multiplicities_index) {

        assert(kit.second.multiplicities.size() == collapsed_best_paths.size());

        uchar max_multiplicity = 0;
        bool has_constant_multiplicity = true;

        for (ushort path_idx = 0; path_idx < kit.second.multiplicities.size(); path_idx++) {

            max_multiplicity = max(max_multiplicity, kit.second.multiplicities.at(path_idx));

            if (kit.second.multiplicities.at(path_idx) != kit.second.multiplicities.front()) {

                has_constant_multiplicity = false;
            }
        }

        assert(max_multiplicity > 0);

        if (max_multiplicity > 1) {

            has_non_unique_kmer = true;
        }
        
        if (kit.second.kmer_counts) {

            auto hash_lock = kmer_hash_hybrid->_hash->lockKey(kit.first);

            static_cast<SampleKmerCounts *>(kit.second.kmer_counts)->addClusterMultiplicity(max_multiplicity, has_constant_multiplicity, variant_cluster_group_idx);

            if (kit.second.kmer_counts->hasMulticlusterOccurrence()) {

                has_non_unique_kmer = true;
      
            } else if (kit.second.kmer_counts->getInterclusterMultiplicity(Utils::Sex::Male) > 0) {

                has_non_unique_kmer = true;
                num_counted_sample_kmers++;  

            } else {

                num_counted_sample_kmers++;                
            }
        }
    }

    for (auto &kmer_multiplicities: path_kmer_multiplicities) {

        if (visited_kmer_vertices.count(kmer_multiplicities) < 1) {

            for (auto &kit: *kmer_multiplicities) {

                if (kit.second.kmer_counts) {

                    auto hash_lock = kmer_hash_hybrid->_hash->lockKey(kit.first);
                    static_cast<SampleKmerCounts *>(kit.second.kmer_counts)->addClusterMultiplicity(0, false, variant_cluster_group_idx);
                }
            }  
        } 

        delete kmer_multiplicities; 
    } 

    return num_counted_sample_kmers;
}


template <uchar kmer_size>
void VariantClusterGraphLockedKmerSize<kmer_size>::getBestHaplotypeCandidates(KmerHash * kmer_hash, VariantClusterHaplotypes * variant_cluster_haplotypes, const uint prng_seed, const ushort num_haplotype_candidates_per_sample) {

    mt19937 prng = mt19937(prng_seed);
    KmerHashHybrid<kmer_size> * kmer_hash_hybrid = static_cast<KmerHashHybrid<kmer_size> * >(kmer_hash);
    
    vector<KmerMultiplicities *> path_kmer_multiplicities;
    path_kmer_multiplicities.reserve(variant_cluster_info.size() * 2 * num_haplotype_candidates_per_sample * (kmer_hash_hybrid->getNumberOfSamples() + 1));

    vector<vector<VariantClusterGraphFullKmerPath<kmer_size> *> > * best_paths_per_sample; 

    best_paths_per_sample = findBestPathsPerSample<VariantClusterGraphFullKmerPath<kmer_size> >(kmer_hash_hybrid, kmer_hash_hybrid->getNumberOfSamples(), &path_kmer_multiplicities, &prng, num_haplotype_candidates_per_sample);   

    assert(best_paths_per_sample->size() == static_cast<uint>(kmer_hash_hybrid->getNumberOfSamples() + 1));

    vector<VariantClusterGraphFullKmerPath<kmer_size> *> collapsed_best_paths = collapsePaths(best_paths_per_sample, num_haplotype_candidates_per_sample * (kmer_hash_hybrid->getNumberOfSamples() + 1));
    delete best_paths_per_sample;

    assert(!(collapsed_best_paths.empty()));
    assert(variant_cluster_haplotypes->empty());

    variant_cluster_haplotypes->variants.reserve(collapsed_best_paths.size());
    variant_cluster_haplotypes->nested_variant_cluster_indices.reserve(collapsed_best_paths.size());
    variant_cluster_haplotypes->haplotype_variant_kmer_indices.reserve(collapsed_best_paths.size());
    variant_cluster_haplotypes->haplotype_multicluster_kmer_indices.reserve(collapsed_best_paths.size());
    variant_cluster_haplotypes->redundant_multicluster_haplotypes.reserve(collapsed_best_paths.size());

    for (ushort path_idx = 0; path_idx < collapsed_best_paths.size(); path_idx++) {

        variant_cluster_haplotypes->variants.emplace_back(unordered_map<ushort, ushort>());
        variant_cluster_haplotypes->nested_variant_cluster_indices.emplace_back(collapsed_best_paths.at(path_idx)->getNestedVariantClusterIndices());
        variant_cluster_haplotypes->haplotype_variant_kmer_indices.emplace_back(unordered_map<ushort, vector<uint> >());
        
        variant_cluster_haplotypes->haplotype_multicluster_kmer_indices.emplace_back(vector<uint>());
        variant_cluster_haplotypes->haplotype_multicluster_kmer_indices.back().reserve(collapsed_best_paths.at(path_idx)->getNumberOfMultiClusterKmers());

        variant_cluster_haplotypes->redundant_multicluster_haplotypes.emplace_back(collapsed_best_paths.size(), true);

        for (auto & variant: collapsed_best_paths.at(path_idx)->getVariants()) {

            assert(variant_cluster_haplotypes->variants.back().emplace(variant.first).second);
            
            auto haplotype_variant_kmer_indices_emplace = variant_cluster_haplotypes->haplotype_variant_kmer_indices.back().emplace(variant.first.first, vector<uint>());
            assert(haplotype_variant_kmer_indices_emplace.second);

            haplotype_variant_kmer_indices_emplace.first->second.reserve(collapsed_best_paths.at(path_idx)->getNumberOfVariantKmers(variant.first.first));
        }
    }

    unordered_set<KmerMultiplicities *> visited_kmer_vertices;
    visited_kmer_vertices.reserve(path_kmer_multiplicities.size());

    KmerMultiplicitiesIndex kmer_multiplicities_index = createKmerMultiplicitiesIndex(&collapsed_best_paths, &visited_kmer_vertices, false);

    for (auto &kmer_multiplicities: path_kmer_multiplicities) {

        delete kmer_multiplicities;
    }

    variant_cluster_haplotypes->kmer_haplotype_multiplicities = Eigen::MatrixXuchar(kmer_multiplicities_index.size(), collapsed_best_paths.size());    
    variant_cluster_haplotypes->kmers.reserve(kmer_multiplicities_index.size());

    variant_cluster_haplotypes->unique_kmer_indices.reserve(kmer_multiplicities_index.size());
    variant_cluster_haplotypes->multicluster_kmer_indices.reserve(kmer_multiplicities_index.size());

    vector<uchar> last_multicluster_kmer_multiplicities(collapsed_best_paths.size(), 0);

    uint cur_row_idx = 0;

    for (auto &kit: kmer_multiplicities_index) {

        assert(kit.second.multiplicities.size() == collapsed_best_paths.size());

        if (kit.second.has_max_multiplcity) {

            assert(!(kit.second.kmer_counts));
            has_excluded_kmer = true;

            continue;
        }
 
        uchar max_multiplicity = 0;
    
        bool has_constant_multiplicity = true;
        bool has_identical_multicluster_kmer_multiplicities = true;

        for (ushort path_idx = 0; path_idx < kit.second.multiplicities.size(); path_idx++) {

            max_multiplicity = max(max_multiplicity, kit.second.multiplicities.at(path_idx));
            variant_cluster_haplotypes->kmer_haplotype_multiplicities(cur_row_idx, path_idx) = kit.second.multiplicities.at(path_idx);

            if (kit.second.multiplicities.at(path_idx) != kit.second.multiplicities.front()) {

                has_constant_multiplicity = false;
            }

            if (kit.second.multiplicities.at(path_idx) != last_multicluster_kmer_multiplicities.at(path_idx)) {

                has_identical_multicluster_kmer_multiplicities = false;
            }
        }

        assert(max_multiplicity > 0);

        if (has_constant_multiplicity) {

            if (kit.second.kmer_counts) { 

                assert(kit.second.kmer_counts->hasConstantMultiplicity());
                assert(kit.second.kmer_counts->isMulti());

                for (ushort sample_idx = 0; sample_idx < kmer_hash_hybrid->getNumberOfSamples(); sample_idx++) {

                    static_cast<MultiClusterKmerCounts *>(kit.second.kmer_counts)->addMulticlusterMultiplicity(sample_idx, max_multiplicity);
                }
            } 

            continue;
        }

        for (auto & path_variants: kit.second.path_variants) {

            assert(kit.second.multiplicities.at(path_variants.path_idx) > 0);

            for (auto & variant_id: path_variants.variant_ids) {

                auto kmer_indices_it = variant_cluster_haplotypes->haplotype_variant_kmer_indices.at(path_variants.path_idx).find(variant_id);

                if (kmer_indices_it != variant_cluster_haplotypes->haplotype_variant_kmer_indices.at(path_variants.path_idx).end()) {

                    kmer_indices_it->second.push_back(variant_cluster_haplotypes->kmers.size());
                }
            }
        }

        if (kit.second.kmer_counts) { 
            
            assert(!(kit.second.kmer_counts)->isExcluded());

            if ((kit.second.kmer_counts)->isMulti()) {

                has_non_unique_kmer = true;
                variant_cluster_haplotypes->multicluster_kmer_indices.push_back(variant_cluster_haplotypes->kmers.size());

                for (ushort path_idx = 0; path_idx < kit.second.multiplicities.size(); path_idx++) {

                    if (kit.second.multiplicities.at(path_idx) > 0) {

                        variant_cluster_haplotypes->haplotype_multicluster_kmer_indices.at(path_idx).push_back(variant_cluster_haplotypes->kmers.size());
                    }

                    if (!has_identical_multicluster_kmer_multiplicities) {

                        for (ushort path_idx_2 = 0; path_idx_2 < path_idx; path_idx_2++) {

                            if (kit.second.multiplicities.at(path_idx) != kit.second.multiplicities.at(path_idx_2)) {

                                variant_cluster_haplotypes->redundant_multicluster_haplotypes.at(path_idx).at(path_idx_2) = false;
                                variant_cluster_haplotypes->redundant_multicluster_haplotypes.at(path_idx_2).at(path_idx) = false;
                            }
                        }
                    }

                    last_multicluster_kmer_multiplicities.at(path_idx) = kit.second.multiplicities.at(path_idx);
                }

            } else {

                variant_cluster_haplotypes->unique_kmer_indices.push_back(variant_cluster_haplotypes->kmers.size());
            }

            variant_cluster_haplotypes->kmers.emplace_back(kit.second.kmer_counts);
        
        } else {

            variant_cluster_haplotypes->unique_kmer_indices.push_back(variant_cluster_haplotypes->kmers.size());
            variant_cluster_haplotypes->kmers.emplace_back(new EmptyKmerCounts());
        }

        cur_row_idx++;
    }

    if (cur_row_idx < variant_cluster_haplotypes->kmer_haplotype_multiplicities.rows()) {

        Eigen::NoChange_t no_change_t;
        variant_cluster_haplotypes->kmer_haplotype_multiplicities.conservativeResize(cur_row_idx, no_change_t);
    }

    variant_cluster_haplotypes->kmers.shrink_to_fit(); 
    variant_cluster_haplotypes->unique_kmer_indices.shrink_to_fit();
    variant_cluster_haplotypes->multicluster_kmer_indices.shrink_to_fit();

    assert(cur_row_idx == variant_cluster_haplotypes->kmer_haplotype_multiplicities.rows());
    assert(cur_row_idx == variant_cluster_haplotypes->kmers.size());
    assert(cur_row_idx == (variant_cluster_haplotypes->unique_kmer_indices.size() + variant_cluster_haplotypes->multicluster_kmer_indices.size()));

    for (auto & variant_kmer_indices: variant_cluster_haplotypes->haplotype_variant_kmer_indices) {

        for (auto & kmer_indices: variant_kmer_indices) {
        
            kmer_indices.second.shrink_to_fit();
        }
    }

    for (auto multicluster_kmer_indices: variant_cluster_haplotypes->haplotype_multicluster_kmer_indices) {

        multicluster_kmer_indices.shrink_to_fit();
    }
}


template <uchar kmer_size>
template <typename PathType>
vector<vector<PathType *> > * VariantClusterGraphLockedKmerSize<kmer_size>::findBestPathsPerSample(KmerHashHybrid<kmer_size> * kmer_hash_hybrid, const ushort num_samples, vector<KmerMultiplicities *> * path_kmer_multiplicities, mt19937 * prng, const ushort num_haplotype_candidates_per_sample) {

    vector<vector<PathType *> > * best_paths_per_sample = nullptr;

    auto vit = boost::vertices(graph);

    auto eit_in = boost::in_edges(*vit.first, graph);
    auto eit_out = boost::out_edges(*vit.first, graph);

    assert(eit_in.first == eit_in.second);
    assert(eit_out.first != eit_out.second);

    unordered_map<vertex_t, vector<vector<PathType *> > * > vertex_best_paths;
    unordered_map<vertex_t, vector<vertex_t> > visited_vertices;

    bitset<2> nt_bits;

    uint expected_path_length = variant_cluster_info.size() * 2;

    while (vit.first != vit.second) {

        eit_in = boost::in_edges(*vit.first, graph);

        if (eit_in.first == eit_in.second) {

            assert(*vit.first == 0);
            assert(graph[*vit.first]->variant_ids.empty());
            assert(!(graph[*vit.first]->isComplex()));

            best_paths_per_sample = new vector<vector<PathType *> >();

            for (ushort sample_idx = 0; sample_idx < num_samples; sample_idx++) {

                best_paths_per_sample->emplace_back(vector<PathType *>(1, new PathType(expected_path_length, false)));
            }

            best_paths_per_sample->emplace_back(vector<PathType *>(1, new PathType(expected_path_length, true)));

        } else {

            vertex_t cur_source_vertex = boost::source(*eit_in.first, graph);

            assert(*vit.first > cur_source_vertex);
            assert(static_cast<uint>(num_samples + 1) == vertex_best_paths.at(cur_source_vertex)->size());

            best_paths_per_sample = new vector<vector<PathType *> >();

            for (auto &sample_paths: *vertex_best_paths.at(cur_source_vertex)) {

                best_paths_per_sample->emplace_back(vector<PathType *>());
                best_paths_per_sample->back().reserve(sample_paths.size());

                for (auto &path: sample_paths) {

                    best_paths_per_sample->back().emplace_back(new PathType(*path));
                }
            } 

            eit_in.first++;
         
            while (eit_in.first != eit_in.second) {

                cur_source_vertex = boost::source(*eit_in.first, graph);

                assert(*vit.first > cur_source_vertex);
                assert(best_paths_per_sample->size() == vertex_best_paths.at(cur_source_vertex)->size());

                for (ushort sample_idx = 0; sample_idx < best_paths_per_sample->size(); sample_idx++) {

                    best_paths_per_sample->at(sample_idx).reserve(best_paths_per_sample->at(sample_idx).size() + vertex_best_paths.at(cur_source_vertex)->at(sample_idx).size());
                    mergePaths(eit_in.first, *vit.first, &best_paths_per_sample->at(sample_idx), &vertex_best_paths.at(cur_source_vertex)->at(sample_idx));
                }

                eit_in.first++;
            }
        }

        uint num_vertex_nucleotides = graph[*vit.first]->getNumberOfNucleotides();

        assert(best_paths_per_sample->size() == static_cast<uint>(num_samples + 1));

        for (auto &sample_paths: *best_paths_per_sample) {

            assert(!(sample_paths.empty()));
            
            sample_paths.shrink_to_fit();
            shuffle(sample_paths.begin(), sample_paths.end(), *prng);

            for (auto &path: sample_paths) {

                path->addVertex(*graph[*vit.first], num_vertex_nucleotides);
            }                 
        }

        bool unique_kmer_vertex = true;
        bool shared_kmer_vertex = false;

        for (auto &sequence: graph[*vit.first]->sequences) {

            if (unique_kmer_vertex) {

                for (auto &sample_paths: *best_paths_per_sample) {

                    for (auto &path: sample_paths) {

                        path_kmer_multiplicities->push_back(path->newKmerVertex(num_vertex_nucleotides));
                    }
                }

                shared_kmer_vertex = true;

            } else if (shared_kmer_vertex) {

                bool is_first = true;

                for (auto &sample_paths: *best_paths_per_sample) {

                    for (auto &path: sample_paths) {

                        if (is_first) {

                            is_first = false;

                            path_kmer_multiplicities->push_back(path->newKmerVertex(num_vertex_nucleotides));
                            path->kmer_pair.reset();

                        } else {

                            path->addKmerVertex(path_kmer_multiplicities->back());
                            path->kmer_pair.reset();
                        }
                    }
                }

                shared_kmer_vertex = false;
            
            } else {

                for (auto &sample_paths: *best_paths_per_sample) {

                    for (auto &path: sample_paths) {

                        path->kmer_pair.reset();
                    }
                }
            }

            auto sit = sequence.begin();

            while (sit != sequence.end()) {

                nt_bits.set(0, *sit);
                sit++;

                nt_bits.set(1, *sit);
                sit++;

                if (unique_kmer_vertex) {

                    for (ushort sample_idx = 0; sample_idx < best_paths_per_sample->size(); sample_idx++) {

                        for (auto &path: best_paths_per_sample->at(sample_idx)) {

                            if (path->kmer_pair.shift(nt_bits)) {    

                                addKmer(path, kmer_hash_hybrid, sample_idx, false);
                            }
                    
                            path->incrementNucleotide();
                        }
                    }

                    if ((((sit - sequence.begin())/2) == (kmer_size - 1)) and (sit != sequence.end())) {

                        bool is_first = true;

                        for (auto &sample_paths: *best_paths_per_sample) {

                            for (auto &path: sample_paths) {

                                if (is_first) {

                                    is_first = false;                                
                                    path_kmer_multiplicities->push_back(path->newKmerVertex(num_vertex_nucleotides));

                                } else {

                                    path->addKmerVertex(path_kmer_multiplicities->back());
                                }
                            }
                        }

                        unique_kmer_vertex = false;
                        shared_kmer_vertex = false;
                    } 

                } else {

                    bool is_first = true;

                    for (ushort sample_idx = 0; sample_idx < best_paths_per_sample->size(); sample_idx++) {

                        for (auto &path: best_paths_per_sample->at(sample_idx)) {

                            if (is_first) {

                                is_first = false;
                            
                                if (path->kmer_pair.shift(nt_bits)) {    

                                    addKmer(path, kmer_hash_hybrid, sample_idx, false);
                                }

                            } else {

                                if (path->kmer_pair.shift(nt_bits)) {    

                                    addKmer(path, kmer_hash_hybrid, sample_idx, true);
                                }
                            }
    
                            path->incrementNucleotide();
                        }
                    }
                }
            }

            unique_kmer_vertex = false;
        }

        assert(best_paths_per_sample->size() == static_cast<uint>(num_samples + 1));     
        
        for (ushort sample_idx = 0; sample_idx < num_samples; sample_idx++) {

            unordered_set<const VariantClusterGraphVertex *> path_covered_vertices;

            assert(!(best_paths_per_sample->at(sample_idx).front()->isDummy()));                 
            filterBestPaths(&(best_paths_per_sample->at(sample_idx)), &path_covered_vertices, static_cast<uint>(num_haplotype_candidates_per_sample * haplotype_candidate_nbest_multiplicity), true);
            assert(!(best_paths_per_sample->at(sample_idx).empty()));            
        }

        unordered_set<const VariantClusterGraphVertex *> path_covered_vertices;

        assert(best_paths_per_sample->back().front()->isDummy());                 
        filterBestPaths(&(best_paths_per_sample->back()), &path_covered_vertices, num_dummy_haplotype_candidates, true);
        assert(!(best_paths_per_sample->back().empty()));            

        assert(vertex_best_paths.insert({*vit.first, best_paths_per_sample}).second);
        deleteVisitedVertexPaths(*vit.first, &vertex_best_paths, &visited_vertices);

        vit.first++;
    }

    assert(vertex_best_paths.size() == 1);
    assert(vertex_best_paths.begin()->second == best_paths_per_sample);

    assert(best_paths_per_sample->size() == static_cast<uint>(num_samples + 1));

    unordered_set<const VariantClusterGraphVertex *> covered_vertices;

    for (ushort sample_idx = 0; sample_idx < num_samples; sample_idx++) {

        unordered_set<const VariantClusterGraphVertex *> path_covered_vertices;

        assert(!(best_paths_per_sample->at(sample_idx).front()->isDummy()));                 
        filterBestPaths(&(best_paths_per_sample->at(sample_idx)), &path_covered_vertices, num_haplotype_candidates_per_sample, true);
        assert(!(best_paths_per_sample->at(sample_idx).empty()));            

        covered_vertices.insert(path_covered_vertices.begin(), path_covered_vertices.end());
    }

    assert(best_paths_per_sample->back().front()->isDummy());                 
    filterBestPaths(&(best_paths_per_sample->back()), &covered_vertices, num_dummy_haplotype_candidates, false);

    return best_paths_per_sample;
}


template <uchar kmer_size>
template <typename PathType>
void VariantClusterGraphLockedKmerSize<kmer_size>::mergePaths(const in_edge_it cur_edge_it, const vertex_t cur_vertex, vector<PathType *> * main_paths, vector<PathType *> * input_paths) {

    assert(cur_vertex > boost::source(*cur_edge_it, graph));
    assert(!(main_paths->empty()));

    auto end_best_path = main_paths->end();
 
    if (hasEqualSourceSequence(boost::in_edges(cur_vertex, graph).first, cur_edge_it)) {

        for (auto &vpit: *input_paths) {

            bool has_equal = false;
            auto bvit = main_paths->begin();

            while (bvit != end_best_path) {

                if (vpit->kmer_pair == (*bvit)->kmer_pair) {

                    if (*vpit == **bvit) {

                        if (swapRedundantPath(**bvit, *vpit)) {

                            delete *bvit;
                            *bvit = new PathType(*vpit);
                        }
                        
                        has_redundant_sequence = true;
                        has_equal = true;
                        break;
                    }  
                } 

                bvit++;                         
            }

            if (!has_equal) {

                main_paths->emplace_back(new PathType(*vpit));
            }
        } 
    
    } else {

        for (auto &vpit: *input_paths) {

            main_paths->emplace_back(new PathType(*vpit));
        }   
    }
}


template <uchar kmer_size>
bool VariantClusterGraphLockedKmerSize<kmer_size>::hasEqualSourceSequence(in_edge_it first_edge, const in_edge_it cur_edge) {

    assert(first_edge != cur_edge);
    auto cur_source_vertex = boost::source(*cur_edge, graph);

    bool has_equal_source_sequence = false;

    while (first_edge != cur_edge) {

        auto fsit = graph[boost::source(*first_edge, graph)]->sequences.back().crbegin();
        auto csit = graph[cur_source_vertex]->sequences.back().crbegin();

        bool has_equal_sequence = true;

        while ((fsit != graph[boost::source(*first_edge, graph)]->sequences.back().crend()) and (csit != graph[cur_source_vertex]->sequences.back().crend())) {

            if (*fsit != *csit) {

                has_equal_sequence = false;
                break;
            }

            fsit++;
            csit++;
        }

        if (has_equal_sequence) {

            has_equal_source_sequence = true;
            break;
        }

        first_edge++;
    }

    return has_equal_source_sequence;
}


template <uchar kmer_size>
template <typename PathType>
bool VariantClusterGraphLockedKmerSize<kmer_size>::swapRedundantPath(PathType & first_redundant_path, PathType & second_redundant_path) {

    if (first_redundant_path.getVariants().size() > second_redundant_path.getVariants().size()) {

        return false;
    
    } else if (second_redundant_path.getVariants().size() > first_redundant_path.getVariants().size()) {

        return true;
    }

    auto first_path_it = first_redundant_path.getPath().begin();
    auto second_path_it = second_redundant_path.getPath().begin();

    while ((first_path_it != first_redundant_path.getPath().end()) and (second_path_it != second_redundant_path.getPath().end())) {

        if (*first_path_it != *second_path_it) {

            assert(!((*first_path_it)->variant_ids.empty()));
            assert(!((*second_path_it)->variant_ids.empty()));
            assert((*first_path_it)->variant_ids.front().is_allele_vertex);
            assert((*second_path_it)->variant_ids.front().is_allele_vertex);
            assert((*first_path_it)->variant_ids.front().allele_id.first == (*second_path_it)->variant_ids.front().allele_id.first);
            assert((*first_path_it)->variant_ids.front().allele_id.second != (*second_path_it)->variant_ids.front().allele_id.second);

            if ((*first_path_it)->variant_ids.front().allele_id.second > (*second_path_it)->variant_ids.front().allele_id.second) {

                return false;
            
            } else {

                return true;
            }
        } 

        first_path_it++;
        second_path_it++;
    }

    assert(false);
    return false;
}


template <uchar kmer_size>
template <typename PathType>
void VariantClusterGraphLockedKmerSize<kmer_size>::addKmer(PathType * path, KmerHashHybrid<kmer_size> * kmer_hash_hybrid, const ushort sample_idx, const bool is_shared) {

    auto skcit = kmer_hash_hybrid->_hash->find(path->kmer_pair.getLexicographicalLowestKmer());

    if (skcit != kmer_hash_hybrid->_hash->end()) {

        assert(*skcit);

        if ((*skcit)->isExcluded()) {

            has_excluded_kmer = true;
        }

        path->newKmer(*skcit, sample_idx, is_shared);

    } else {

        path->newKmer(nullptr, sample_idx, is_shared);
    }
}


template <uchar kmer_size>
template <typename PathType>
void VariantClusterGraphLockedKmerSize<kmer_size>::filterBestPaths(vector<PathType *> * best_paths, unordered_set<const VariantClusterGraphVertex *> * path_covered_vertices, const uint max_paths, const bool keep_zero_covered_paths) {

    assert(!(best_paths->empty()));

    if ((best_paths->size() > max_paths) or (!keep_zero_covered_paths)) {

        assert(best_paths->size() > 1);

        vector<PathType *> positive_cover_sorted_paths;
        vector<PathType *> zero_cover_sorted_paths;

        unordered_set<PathType *> scored_paths;

        const uint paths_to_score = min(static_cast<uint>(best_paths->size()), max_paths);
        assert(paths_to_score > 0);

        while (scored_paths.size() < paths_to_score) {

            pair<pair<double, uint>, PathType *> highest_scoring_path = make_pair(make_pair(-1, 0), nullptr);

            for (auto &pit: *best_paths) {

                if (scored_paths.count(pit) < 1) {

                    pair<double, uint> cur_path_score = pit->getScore(*path_covered_vertices);

                    if (Utils::doubleCompare(cur_path_score.first, highest_scoring_path.first.first)) {

                        if (cur_path_score.second > highest_scoring_path.first.second) {

                            highest_scoring_path.first = cur_path_score;
                            highest_scoring_path.second = pit;
                        }
                   
                    } else if (cur_path_score.first > highest_scoring_path.first.first) {

                        highest_scoring_path.first = cur_path_score;
                        highest_scoring_path.second = pit;                   
                    }
                }
            }

            assert(highest_scoring_path.first.first >= 0);
            assert(highest_scoring_path.second);

            if (highest_scoring_path.first.second > 0) {

                highest_scoring_path.second->addCoveredVertices(path_covered_vertices);
                positive_cover_sorted_paths.push_back(highest_scoring_path.second);    

            } else {

                zero_cover_sorted_paths.push_back(highest_scoring_path.second);
            }

            assert(scored_paths.insert(highest_scoring_path.second).second);
        }

        assert(path_covered_vertices->size() <= boost::num_vertices(graph));
        
        assert(scored_paths.size() == paths_to_score);
        assert(scored_paths.size() == (positive_cover_sorted_paths.size() + zero_cover_sorted_paths.size()));

        const bool is_dummy = best_paths->front()->isDummy();

        if (!keep_zero_covered_paths) {

            assert(is_dummy);
        }

        for (auto &path: *best_paths) {
            
            assert(path->isDummy() == is_dummy);

            if (scored_paths.count(path) < 1) {

                delete path;
            }
        }

        best_paths->clear();
        
        for (auto &path: positive_cover_sorted_paths) {

            best_paths->push_back(path);
        }

        if (keep_zero_covered_paths) {

            for (auto &path: zero_cover_sorted_paths) {

                best_paths->push_back(path);
            }           
        
        } else {

           for (auto &path: zero_cover_sorted_paths) {

                delete path;
            }   
        }
    }

    assert(best_paths->size() <= max_paths);
}


template <uchar kmer_size>
template <typename PathType>
void VariantClusterGraphLockedKmerSize<kmer_size>::deleteVisitedVertexPaths(vertex_t cur_vertex, unordered_map<vertex_t, vector<vector<PathType *> > * > * vertex_best_paths, unordered_map<vertex_t, vector<vertex_t> > * visited_vertices) {
 
    auto eit_out = boost::out_edges(cur_vertex, graph); 
    vertex_t max_vertex_id = 0;

    while (eit_out.first != eit_out.second) {

        vertex_t vertex_id = boost::target(*eit_out.first, graph);
        assert(vertex_id > cur_vertex);

        max_vertex_id = max(max_vertex_id, vertex_id);
        eit_out.first++;
    }

    if (!(visited_vertices->insert({max_vertex_id, vector<vertex_t>(1, cur_vertex)}).second)) {

        visited_vertices->at(max_vertex_id).push_back(cur_vertex);
    }

    auto visited_vertices_it = visited_vertices->find(cur_vertex);

    if (visited_vertices_it != visited_vertices->end()) {

        for (auto lit: visited_vertices_it->second) {

            auto vertex_best_paths_it = vertex_best_paths->find(lit);
            assert(vertex_best_paths_it != vertex_best_paths->end());

            for (auto &sample_paths: *vertex_best_paths_it->second) {

                for (auto &path: sample_paths) {

                    delete path;
                }
            }

            delete vertex_best_paths_it->second;
            assert(vertex_best_paths->erase(lit));
        }           
    }
}
      

template <uchar kmer_size>
template <typename PathType>
vector<PathType *> VariantClusterGraphLockedKmerSize<kmer_size>::collapsePaths(vector<vector<PathType *> > * best_paths_per_sample, const uint max_collapsed_paths) {

    vector<PathType *> collapsed_best_paths;
    collapsed_best_paths.reserve(max_collapsed_paths);

    for (auto &sample_paths: *best_paths_per_sample) {

        auto sample_paths_it = sample_paths.begin();

        while (sample_paths_it != sample_paths.end()) {

            assert(*sample_paths_it);
            bool has_equal = false;

            auto collapsed_best_paths_it = collapsed_best_paths.begin();

            while (collapsed_best_paths_it != collapsed_best_paths.end()) {

                if ((*sample_paths_it)->kmer_pair == (*collapsed_best_paths_it)->kmer_pair) {

                    if ((*sample_paths_it)->getPath() == (*collapsed_best_paths_it)->getPath()) {

                        has_equal = true;
                        break;

                    } else if (**sample_paths_it == **collapsed_best_paths_it) {

                        if (swapRedundantPath(**collapsed_best_paths_it, **sample_paths_it)) {

                            auto temp_path_ptr = *collapsed_best_paths_it;
                            
                            *collapsed_best_paths_it = *sample_paths_it;
                            *sample_paths_it = temp_path_ptr;
                        }

                        has_equal = true;
                        has_redundant_sequence = true;
                        break;
                    }
                }

                collapsed_best_paths_it++;
            }

            if (has_equal) {

                delete *sample_paths_it;
            
            } else {

                collapsed_best_paths.push_back(*sample_paths_it);
            }

            sample_paths_it++;
        }
    } 

    collapsed_best_paths.shrink_to_fit();
    return collapsed_best_paths; 
} 


template<uchar kmer_size>
template <typename PathType>
typename VariantClusterGraphLockedKmerSize<kmer_size>::KmerMultiplicitiesIndex VariantClusterGraphLockedKmerSize<kmer_size>::createKmerMultiplicitiesIndex(vector<PathType *> * collapsed_best_paths, unordered_set<KmerMultiplicities *> * visited_kmer_vertices, const bool is_pre_counting) {

    KmerMultiplicitiesIndex kmer_multiplicities_index;

    for (ushort path_idx = 0; path_idx < collapsed_best_paths->size(); path_idx++) {

        for (auto &kmer_vertex: collapsed_best_paths->at(path_idx)->getKmerVertexMultiplicities()) {

            visited_kmer_vertices->insert(kmer_vertex);

            for (auto &kmer: *kmer_vertex) {

                auto kmer_multiplicities_index_emplace = kmer_multiplicities_index.emplace(kmer.first, PathsKmerInfo<kmer_size>(kmer.second.kmer_counts, collapsed_best_paths->size()));

                assert(kmer_multiplicities_index_emplace.first->second.kmer_counts == kmer.second.kmer_counts);

                if (kmer_multiplicities_index_emplace.first->second.multiplicities.at(path_idx) < Utils::uchar_overflow) {

                    kmer_multiplicities_index_emplace.first->second.multiplicities.at(path_idx)++;
                
                } else {

                    kmer_multiplicities_index_emplace.first->second.has_max_multiplcity = true;
                }

                if (!is_pre_counting) {

                    vector<PathVariants> * path_variants = &(kmer_multiplicities_index_emplace.first->second.path_variants);

                    if (path_variants->empty()) {

                        path_variants->emplace_back(path_idx);
                        path_variants->back().variant_ids = kmer.second.variant_ids;
                    
                    } else if (path_variants->back().path_idx != path_idx) {

                        path_variants->emplace_back(path_idx);
                        path_variants->back().variant_ids = kmer.second.variant_ids;                      
                    
                    } else {

                        for (auto & variant_id: kmer.second.variant_ids) {

                            if (find(path_variants->back().variant_ids.begin(), path_variants->back().variant_ids.end(), variant_id) == path_variants->back().variant_ids.end()) {

                                path_variants->back().variant_ids.push_back(variant_id);
                            }  
                        }
                    } 
                }
            }       
        }

        delete collapsed_best_paths->at(path_idx);
    }

    return kmer_multiplicities_index;
}



template class VariantClusterGraphLockedKmerSize<31>;
template class VariantClusterGraphLockedKmerSize<39>;
template class VariantClusterGraphLockedKmerSize<47>;
template class VariantClusterGraphLockedKmerSize<55>;
template class VariantClusterGraphLockedKmerSize<63>;

