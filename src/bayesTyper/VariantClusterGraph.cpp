
/*
VariantClusterGraph.cpp - This file is part of BayesTyper (v1.1)


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
#include <queue>

#include "boost/graph/adjacency_list.hpp"

#include "VariantClusterGraph.hpp"
#include "Utils.hpp"
#include "VariantCluster.hpp"
#include "KmerHash.hpp"
#include "PerfectSet.hpp"
#include "VariantClusterGraphPath.hpp"
#include "Kmer.hpp"
#include "Sequence.hpp"
#include "KmerCounts.hpp"
#include "VariantClusterHaplotypes.hpp"
#include "VariantInfo.hpp"
#include "VariantClusterGraphVertex.hpp"

static const ushort max_smallmer_counter_recursion_depth = 100;
static const uint max_total_paths = 256;


VariantClusterGraph::VariantClusterGraph(VariantCluster * variant_cluster, const string & chromosome_sequence, const uchar kmer_size) {

    has_ambiguous_nucleotide = false;
    has_redundant_sequence = false;
    has_intercluster_kmer = false;
    has_excluded_kmer = false;

    assert(variant_cluster->variants.size() < Utils::ushort_overflow);
    
    assert(variant_cluster_info.empty());
    variant_cluster_info.reserve(variant_cluster->variants.size());

    map<uint, pair<vector<vertex_t>, vector<ushort> > > added_vertices;
    unordered_set<ushort> reference_variant_indices;

    auto lit = variant_cluster->variants.begin();
    auto cit = chromosome_sequence.begin();

    vertex_t cur_vertex = boost::add_vertex(graph);    

    addVertices(&cur_vertex, vector<StringItPair>(1, StringItPair(cit + lit->first - (kmer_size - 1), cit + lit->first)), make_pair(Utils::ushort_overflow, Utils::ushort_overflow), reference_variant_indices, vector<uint>());

    vertex_t prev_vertex = cur_vertex;

    assert(added_vertices.insert({lit->first, make_pair(vector<vertex_t>(1, cur_vertex), vector<ushort>())}).second);
    pair<uint, vertex_t> prev_new_position(lit->first, cur_vertex);

    uint cur_last_position = 0;
    uint next_position = 0;

    ushort variant_counter = 0;

    while (lit != variant_cluster->variants.end()) {

        variant_cluster_info.emplace_back(VariantInfo());
        variant_cluster_info.back().variant_id = lit->second.id;
        variant_cluster_info.back().num_alleles = lit->second.alternative_alleles.size() + 1 + lit->second.excluded_alt_alleles.size() + lit->second.has_dependency;
        variant_cluster_info.back().excluded_alt_alleles = lit->second.excluded_alt_alleles;
        variant_cluster_info.back().has_dependency = lit->second.has_dependency;

        for (auto & ait: lit->second.alternative_alleles) {

            assert(count(lit->second.excluded_alt_alleles.begin(), lit->second.excluded_alt_alleles.end(), ait.first) == 0);

            vertex_t next_vertex = boost::add_vertex(graph);

            if (prev_new_position.first == lit->first) { 

                assert(boost::add_edge(prev_new_position.second, next_vertex, graph).second);

            } else {

                assert(boost::add_edge(cur_vertex, next_vertex, graph).second);
            }

            addVertices(&next_vertex, vector<StringItPair>(1, StringItPair(ait.second.second.begin(), ait.second.second.end())), make_pair(variant_counter, ait.first), reference_variant_indices, vector<uint>());

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

        reference_variant_indices.insert(variant_counter);

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
                
                reference_variant_indices.erase(variant_id);
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

                addVertices(&cur_vertex, contained_vertices, pair<ushort, ushort>(variant_counter, 0), reference_variant_indices, nested_variant_cluster_indices);

            } else {

                addVertices(&cur_vertex, contained_vertices, pair<ushort, ushort>(Utils::ushort_overflow, Utils::ushort_overflow), reference_variant_indices, nested_variant_cluster_indices);
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

void VariantClusterGraph::addVertices(vertex_t * cur_vertex, const vector<StringItPair> & vertex_sequences, const pair<ushort, ushort> & variant_allele_idx, const unordered_set<ushort> & reference_variant_indices, const vector<uint> & nested_variant_cluster_indices) {

    assert(!(vertex_sequences.empty()));
    assert(vertex_sequences.size() == (nested_variant_cluster_indices.size() + 1));

    vector<ushort> vertex_reference_variant_indices;
    vertex_reference_variant_indices.reserve(reference_variant_indices.size());

    for (auto & reference_variant_idx: reference_variant_indices) {

        assert(reference_variant_idx != Utils::ushort_overflow);
    
        if (reference_variant_idx != variant_allele_idx.first) {

            vertex_reference_variant_indices.push_back(reference_variant_idx);
        }
    }

    initVertex(cur_vertex, vertex_sequences.front(), variant_allele_idx, vertex_reference_variant_indices, Utils::uint_overflow);
    
    for (uint vertex_idx = 1; vertex_idx < vertex_sequences.size(); vertex_idx++) {

        vertex_t prev_vertex = *cur_vertex;
        *cur_vertex = boost::add_vertex(graph);

        assert(boost::add_edge(prev_vertex, *cur_vertex, graph).second);

        initVertex(cur_vertex, vertex_sequences.at(vertex_idx), variant_allele_idx, vertex_reference_variant_indices, nested_variant_cluster_indices.at(vertex_idx - 1));
    }
}

void VariantClusterGraph::initVertex(vertex_t * cur_vertex, StringItPair vertex_sequence, const pair<ushort, ushort> & variant_allele_idx, const vector<ushort> & vertex_reference_variant_indices, const uint nested_variant_cluster_index) {

    graph[*cur_vertex].variant_allele_idx = variant_allele_idx;
    graph[*cur_vertex].reference_variant_indices = vertex_reference_variant_indices;
    graph[*cur_vertex].nested_variant_cluster_index = nested_variant_cluster_index;

    if (nested_variant_cluster_index != Utils::uint_overflow) {

        graph[*cur_vertex].is_disconnected = true;

    } else {

        graph[*cur_vertex].is_disconnected = false;
    }

    graph[*cur_vertex].sequence.reserve((vertex_sequence.second - vertex_sequence.first) * 2);

    bool prev_is_disconnected = false;

    while (vertex_sequence.first != vertex_sequence.second) {

        auto nt_bits = Sequence::ntToBit<1>(*vertex_sequence.first);

        if (!nt_bits.second) {

            has_ambiguous_nucleotide = true;

            if (!prev_is_disconnected) {

                graph[*cur_vertex].sequence.shrink_to_fit();
                assert((graph[*cur_vertex].sequence.size() % 2) == 0);

                vertex_t prev_vertex = *cur_vertex;
                *cur_vertex = boost::add_vertex(graph);

                assert(boost::add_edge(prev_vertex, *cur_vertex, graph).second);

                graph[*cur_vertex].variant_allele_idx = variant_allele_idx;
                graph[*cur_vertex].reference_variant_indices = vertex_reference_variant_indices;
                graph[*cur_vertex].nested_variant_cluster_index = nested_variant_cluster_index;
                graph[*cur_vertex].is_disconnected = true;
                graph[*cur_vertex].sequence.reserve((vertex_sequence.second - vertex_sequence.first) * 2);
            }

            prev_is_disconnected = true;

        } else {

            graph[*cur_vertex].sequence.push_back(nt_bits.first[0]);
            graph[*cur_vertex].sequence.push_back(nt_bits.first[1]);

            prev_is_disconnected = false;
        }
        
        vertex_sequence.first++;
    }

    assert((graph[*cur_vertex].sequence.size() % 2) == 0);
}

ulong VariantClusterGraph::countSmallmers(Utils::SmallmerSet * smallmer_set) {

    ulong num_unique_smallmers = 0;

    auto vit = boost::vertices(graph);

    auto eit_in = boost::in_edges(*vit.first, graph);
    auto eit_out = boost::out_edges(*vit.first, graph);

    assert(eit_in.first == eit_in.second);
    assert(eit_out.first != eit_out.second);

    set<vertex_t> remaining_vertices;

    while (vit.first != vit.second) {

        assert(remaining_vertices.insert(*vit.first).second);
        vit.first++;
    }

    assert(!(remaining_vertices.empty()));

    auto remaining_vertices_it = remaining_vertices.begin();

    while (remaining_vertices_it != remaining_vertices.end()) {

        auto vertex = *remaining_vertices_it;
        remaining_vertices.erase(remaining_vertices_it);

        KmerPair<Utils::small_kmer_size> kmer_pair;
        num_unique_smallmers += countVertexSmallmers(vertex, smallmer_set, kmer_pair, &remaining_vertices, Utils::uint_overflow, 0);

        remaining_vertices_it = remaining_vertices.begin();
    }

    return num_unique_smallmers;
}


ulong VariantClusterGraph::countVertexSmallmers(const vertex_t vertex, Utils::SmallmerSet * smallmer_set, KmerPair<Utils::small_kmer_size> kmer_pair, set<vertex_t> * remaining_vertices, uint remaining_prefix_size, ushort recursion_depth) {

    ulong num_unique_smallmers = 0;
    recursion_depth++;

    assert((remaining_prefix_size < Utils::small_kmer_size) or (remaining_prefix_size == Utils::uint_overflow));

    if (graph[vertex].is_disconnected) {

        kmer_pair.reset();
    }

    bitset<2> nt_bits;

    assert((graph[vertex].sequence.size() % 2) == 0);    
    auto sequence_it = graph[vertex].sequence.cbegin();

    while (sequence_it != graph[vertex].sequence.cend()) {

        nt_bits.set(0, *sequence_it);
        sequence_it++;

        nt_bits.set(1, *sequence_it);
        sequence_it++;

       if (kmer_pair.move(nt_bits)) { 

            num_unique_smallmers += smallmer_set->insert(kmer_pair.getForwardKmer());
            num_unique_smallmers += smallmer_set->insert(kmer_pair.getReverseComplementKmer());
        } 

        remaining_prefix_size--;

        if (remaining_prefix_size == 0) {

            break;
        } 
    }

    if (remaining_prefix_size > 0) {

        auto out_eit = boost::out_edges(vertex, graph);

        while (out_eit.first != out_eit.second) {

            auto target_vertex = boost::target(*out_eit.first, graph);
            auto remaining_vertices_it = remaining_vertices->find(target_vertex);
                        
            if (remaining_prefix_size < Utils::small_kmer_size) {

                num_unique_smallmers += countVertexSmallmers(target_vertex, smallmer_set, kmer_pair, remaining_vertices, remaining_prefix_size, recursion_depth);

            } else if (remaining_vertices_it == remaining_vertices->end()) {

                num_unique_smallmers += countVertexSmallmers(target_vertex, smallmer_set, kmer_pair, remaining_vertices, Utils::small_kmer_size - 1, recursion_depth);

            } else if (recursion_depth >= max_smallmer_counter_recursion_depth) {

                num_unique_smallmers += countVertexSmallmers(target_vertex, smallmer_set, kmer_pair, remaining_vertices, Utils::small_kmer_size - 1, recursion_depth);
            
            } else {

                remaining_vertices->erase(remaining_vertices_it);
                num_unique_smallmers += countVertexSmallmers(target_vertex, smallmer_set, kmer_pair, remaining_vertices, Utils::uint_overflow, recursion_depth);
            } 

            out_eit.first++;
        }   
    }     

    return num_unique_smallmers;
}

const vector<VariantInfo> & VariantClusterGraph::getInfo() {

    return variant_cluster_info;  
}

bool VariantClusterGraph::hasAmbiguousNucleotide() {

    return has_ambiguous_nucleotide;
}

bool VariantClusterGraph::hasRedundantSequence() {

    return has_redundant_sequence;
}

bool VariantClusterGraph::hasInterclusterKmer() {

    return has_intercluster_kmer;
}

bool VariantClusterGraph::hasExcludedKmer() {

    return has_excluded_kmer;
}


template <uchar kmer_size>
VariantClusterKmerGraph<kmer_size>::VariantClusterKmerGraph(VariantCluster * variant_cluster, const string & chromosome_sequence) : VariantClusterGraph(variant_cluster, chromosome_sequence, kmer_size) {
}

template <uchar kmer_size>
void VariantClusterKmerGraph<kmer_size>::countKmers(KmerHash * kmer_hash, const uint variant_cluster_group_idx, const uint prng_seed, const ushort num_samples, const ushort max_sample_haplotype_candidates) {
    
    mt19937 prng = mt19937(prng_seed);

    assert(best_paths.empty());
    findBestPaths(kmer_hash, &prng, num_samples, max_sample_haplotype_candidates);   

    const uint num_vertices = boost::num_vertices(graph);

    assert(!(best_paths.empty()));
    assert((best_paths.size() % num_vertices) == 0);
    assert(best_paths.size() <= (num_vertices * max_total_paths));

    const ushort num_best_paths = best_paths.size() / num_vertices;

    auto kmer_multiplicities_index = indexKmerMultiplicities(kmer_hash);

    for (auto &kit: kmer_multiplicities_index.index) {

        assert(kit.second.kmer_counts);
        assert(kit.second.multiplicities.size() == num_best_paths);

        if ((kit.second.kmer_counts)->getInterclusterMultiplicity(Utils::Gender::Male) > 0) {

            has_intercluster_kmer = true;
        }

        uchar max_multiplicity = 0;
        bool has_constant_multiplicity = true;

        for (ushort path_idx = 0; path_idx < kit.second.multiplicities.size(); path_idx++) {

            max_multiplicity = max(max_multiplicity, kit.second.multiplicities.at(path_idx));

            if (kit.second.multiplicities.at(path_idx) != kit.second.multiplicities.front()) {

                has_constant_multiplicity = false;
            }
        }

        assert(max_multiplicity > 0);
        
        auto hash_lock = static_cast<BasicKmerHash<kmer_size> *>(kmer_hash)->getKmerLock(kit.first);
        kit.second.kmer_counts->addClusterMultiplicity(max_multiplicity, has_constant_multiplicity, variant_cluster_group_idx);
    }
}

template <uchar kmer_size>
void VariantClusterKmerGraph<kmer_size>::getBestHaplotypeCandidates(KmerHash * kmer_hash, VariantClusterHaplotypes * variant_cluster_haplotypes, const ushort num_samples, const ushort max_sample_haplotype_candidates, const uchar num_genomic_rate_gc_bias_bins) {

    const uint num_vertices = boost::num_vertices(graph);

    assert(!(best_paths.empty()));
    assert((best_paths.size() % num_vertices) == 0);
    
    const ushort num_best_paths = best_paths.size() / num_vertices;

    assert(variant_cluster_haplotypes->empty());
    variant_cluster_haplotypes->haplotypes.reserve(best_paths.size());

    for (ushort path_idx = 0; path_idx < num_best_paths; path_idx++) {

        variant_cluster_haplotypes->haplotypes.emplace_back(variant_cluster_info.size());

        for (uint vertex_idx = 0; vertex_idx < num_vertices; vertex_idx++) {

            if (best_paths.at((path_idx * num_vertices) + vertex_idx)) {

                if (graph[vertex_idx].variant_allele_idx.first != Utils::ushort_overflow) {

                    assert(graph[vertex_idx].variant_allele_idx.second != Utils::ushort_overflow);

                    if (!(graph[vertex_idx].is_disconnected)) {

                        assert(variant_cluster_haplotypes->haplotypes.back().variant_allele_indices.at(graph[vertex_idx].variant_allele_idx.first) == Utils::ushort_overflow);
                        variant_cluster_haplotypes->haplotypes.back().variant_allele_indices.at(graph[vertex_idx].variant_allele_idx.first) = graph[vertex_idx].variant_allele_idx.second;
                    }

                    assert(variant_cluster_haplotypes->haplotypes.back().variant_allele_indices.at(graph[vertex_idx].variant_allele_idx.first) == graph[vertex_idx].variant_allele_idx.second);
                }

                if (graph[vertex_idx].nested_variant_cluster_index != Utils::uint_overflow) {

                    assert(graph[vertex_idx].is_disconnected);
                    variant_cluster_haplotypes->haplotypes.back().nested_variant_cluster_indices.push_back(graph[vertex_idx].nested_variant_cluster_index);
                }
            }
        }

        sort(variant_cluster_haplotypes->haplotypes.back().nested_variant_cluster_indices.begin(), variant_cluster_haplotypes->haplotypes.back().nested_variant_cluster_indices.end());

        for (ushort variant_idx = 0; variant_idx < variant_cluster_info.size(); variant_idx++) {

            if (variant_cluster_haplotypes->haplotypes.back().variant_allele_indices.at(variant_idx) == Utils::ushort_overflow) {

                assert(variant_cluster_info.at(variant_idx).has_dependency);
                variant_cluster_haplotypes->haplotypes.back().variant_allele_indices.at(variant_idx) = variant_cluster_info.at(variant_idx).num_alleles - 1;
            }
        }        
    }

    auto variant_kmer_multiplicities_index = indexVariantKmerMultiplicities(kmer_hash);

    variant_cluster_haplotypes->haplotype_unique_kmer_multiplicities = Eigen::MatrixXuchar(variant_kmer_multiplicities_index.num_unique_kmers, num_best_paths); 
    variant_cluster_haplotypes->haplotype_multicluster_kmer_multiplicities = Eigen::MatrixXuchar(variant_kmer_multiplicities_index.num_multicluster_kmers, num_best_paths); 

    variant_cluster_haplotypes->unique_kmers.reserve(variant_kmer_multiplicities_index.num_unique_kmers);
    variant_cluster_haplotypes->multicluster_kmers.reserve(variant_kmer_multiplicities_index.num_multicluster_kmers);

    uint unique_row_idx = 0;
    uint multicluster_row_idx = 0;

    for (auto &kit: variant_kmer_multiplicities_index.index) {

        assert(kit.second.multiplicities.size() == num_best_paths);
        
        bool is_unique = true;

        if (kit.second.kmer_counts) { 

            assert((kit.second.kmer_counts)->hasClusterOccurrence());
            assert(!(kit.second.kmer_counts)->isExcluded());
        
            if ((kit.second.kmer_counts)->hasMulticlusterOccurrence()) {

                assert(variant_kmer_multiplicities_index.num_multicluster_kmers > 0);
                is_unique = false; 
            }
        }

        uchar max_multiplicity = 0;
        bool has_constant_multiplicity = true;

        for (ushort path_idx = 0; path_idx < kit.second.multiplicities.size(); path_idx++) {

            max_multiplicity = max(max_multiplicity, kit.second.multiplicities.at(path_idx));

            if (is_unique) {

                variant_cluster_haplotypes->haplotype_unique_kmer_multiplicities(unique_row_idx, path_idx) = kit.second.multiplicities.at(path_idx);

            } else {

                variant_cluster_haplotypes->haplotype_multicluster_kmer_multiplicities(multicluster_row_idx, path_idx) = kit.second.multiplicities.at(path_idx);                
            }

            if (kit.second.multiplicities.at(path_idx) != kit.second.multiplicities.front()) {

                has_constant_multiplicity = false;
            }
        }

        assert(max_multiplicity > 0);

        if (max_multiplicity > Utils::bit7_overflow) {

            assert(!(kit.second.kmer_counts));
            
            has_excluded_kmer = true;
            continue;
        }

        if (has_constant_multiplicity) {

            if (kit.second.kmer_counts) {

                assert(!(kit.second.kmer_counts->hasConstantMultiplicity()));
                assert(kit.second.kmer_counts->hasMulticlusterOccurrence());
            
            } else {

                has_excluded_kmer = true;
                continue;
            }
        }

        if (is_unique) { 

            variant_cluster_haplotypes->unique_kmers.emplace_back(kit.second.kmer_counts, Sequence::gcBiasBin<kmer_size>(kit.first, num_genomic_rate_gc_bias_bins), kit.second.variant_path_indices);
            unique_row_idx++;

        } else {

            variant_cluster_haplotypes->multicluster_kmers.emplace_back(kit.second.kmer_counts, Sequence::gcBiasBin<kmer_size>(kit.first, num_genomic_rate_gc_bias_bins), kit.second.variant_path_indices);
            multicluster_row_idx++;
        }
    }

    if (unique_row_idx < variant_cluster_haplotypes->haplotype_unique_kmer_multiplicities.rows()) {

        Eigen::NoChange_t no_change_t;
        variant_cluster_haplotypes->haplotype_unique_kmer_multiplicities.conservativeResize(unique_row_idx, no_change_t);
    }

    if (multicluster_row_idx < variant_cluster_haplotypes->haplotype_multicluster_kmer_multiplicities.rows()) {

        Eigen::NoChange_t no_change_t;
        variant_cluster_haplotypes->haplotype_multicluster_kmer_multiplicities.conservativeResize(multicluster_row_idx, no_change_t);
    }

    variant_cluster_haplotypes->unique_kmers.shrink_to_fit(); 
    variant_cluster_haplotypes->multicluster_kmers.shrink_to_fit();

    assert(unique_row_idx == variant_cluster_haplotypes->haplotype_unique_kmer_multiplicities.rows());
    assert(unique_row_idx == variant_cluster_haplotypes->unique_kmers.size());

    assert(multicluster_row_idx == variant_cluster_haplotypes->haplotype_multicluster_kmer_multiplicities.rows());
    assert(multicluster_row_idx == variant_cluster_haplotypes->multicluster_kmers.size());
}

template <uchar kmer_size>
void VariantClusterKmerGraph<kmer_size>::findBestPaths(KmerHash * kmer_hash, mt19937 * prng, const ushort num_samples, const ushort max_sample_paths) {

    auto vit = boost::vertices(graph);

    auto eit_in = boost::in_edges(*vit.first, graph);
    auto eit_out = boost::out_edges(*vit.first, graph);

    assert(eit_in.first == eit_in.second);
    assert(eit_out.first != eit_out.second);

    unordered_map<vertex_t, vector<vector<VariantClusterGraphPath<kmer_size> *> > > vertex_best_paths;
    unordered_map<vertex_t, vertex_t> visited_vertices;

    vector<vector<VariantClusterGraphPath<kmer_size> *> > * best_paths_per_sample = nullptr;

    while (vit.first != vit.second) {

        auto vertex_best_paths_it = vertex_best_paths.emplace(*vit.first, vector<vector<VariantClusterGraphPath<kmer_size> *> >(num_samples));
        assert(vertex_best_paths_it.second);
    
        best_paths_per_sample = &(vertex_best_paths_it.first->second);
        assert(best_paths_per_sample);

        eit_in = boost::in_edges(*vit.first, graph);

        if (eit_in.first == eit_in.second) {

            assert(*vit.first == 0);

            assert(graph[*vit.first].variant_allele_idx.first == Utils::ushort_overflow);
            assert(graph[*vit.first].variant_allele_idx.second == Utils::ushort_overflow);

            assert(graph[*vit.first].reference_variant_indices.empty());
            assert(graph[*vit.first].nested_variant_cluster_index == Utils::uint_overflow);

            for (ushort sample_idx = 0; sample_idx < num_samples; sample_idx++) {

                best_paths_per_sample->at(sample_idx).emplace_back(new VariantClusterGraphPath<kmer_size>(variant_cluster_info.size()));
            }

        } else {
         
            while (eit_in.first != eit_in.second) {

                auto cur_source_vertex = boost::source(*eit_in.first, graph);

                assert(cur_source_vertex < *vit.first);
                assert(best_paths_per_sample->size() == vertex_best_paths.at(cur_source_vertex).size());

                auto visited_vertices_it = visited_vertices.find(cur_source_vertex);
            
                assert(visited_vertices_it != visited_vertices.end());
                assert(visited_vertices_it->first != visited_vertices_it->second);

                for (ushort sample_idx = 0; sample_idx < best_paths_per_sample->size(); sample_idx++) {

                    mergePaths(&(best_paths_per_sample->at(sample_idx)), &(vertex_best_paths.at(cur_source_vertex).at(sample_idx)), *vit.first == visited_vertices_it->second);
                }

                eit_in.first++;
            }
        }

        assert(best_paths_per_sample);
        assert(best_paths_per_sample->size() == num_samples);

        for (auto & sample_paths: *best_paths_per_sample) {
            
            sample_paths.shrink_to_fit();
            shuffle(sample_paths.begin(), sample_paths.end(), *prng);

            assert(!(sample_paths.empty()));

            for (auto & path: sample_paths) {

                if (graph[*vit.first].is_disconnected) {

                    path->kmer_pair.reset();
                }

                path->addVertex(*vit.first, graph[*vit.first]);
            }
        }

        bitset<2> nt_bits;

        assert((graph[*vit.first].sequence.size() % 2) == 0);
        auto sequence_it = graph[*vit.first].sequence.begin();

        while (sequence_it != graph[*vit.first].sequence.end()) {

            nt_bits.set(0, *sequence_it);
            sequence_it++;

            nt_bits.set(1, *sequence_it);
            sequence_it++;

            for (ushort sample_idx = 0; sample_idx < best_paths_per_sample->size(); sample_idx++) {

                for (auto & path: best_paths_per_sample->at(sample_idx)) {

                    if (path->kmer_pair.move(nt_bits)) {    

                        path->updateScore(static_cast<BasicKmerHash<kmer_size> *>(kmer_hash)->findKmer(path->kmer_pair.getLexicographicalLowestKmer()), sample_idx);
                    }
                }
            }
        }

        assert(best_paths_per_sample->size() == num_samples);     
        
        for (ushort sample_idx = 0; sample_idx < num_samples; sample_idx++) {

            filterBestPaths(&(best_paths_per_sample->at(sample_idx)), max_sample_paths, false);
            assert(!(best_paths_per_sample->at(sample_idx).empty()));            
        }       

        vertex_t max_visited_vertices = *vit.first;
        eit_out = boost::out_edges(*vit.first, graph);

        while (eit_out.first != eit_out.second) {

            vertex_t target_vertex_id = boost::target(*eit_out.first, graph);
            assert(*vit.first < target_vertex_id);

            max_visited_vertices = max(max_visited_vertices, target_vertex_id);
            eit_out.first++;
        }

        assert(visited_vertices.emplace(*vit.first, max_visited_vertices).second);

        vit.first++;
    }

    assert(best_paths_per_sample);
    assert(best_paths_per_sample->size() == num_samples);

    for (ushort sample_idx = 0; sample_idx < num_samples; sample_idx++) {

        filterBestPaths(&(best_paths_per_sample->at(sample_idx)), max_sample_paths, true);
        assert(!(best_paths_per_sample->at(sample_idx).empty()));            
    }

    auto collapsed_best_paths = collapsePaths(best_paths_per_sample, prng, max_sample_paths);

    const uint num_vertices = boost::num_vertices(graph);

    assert(best_paths.empty());
    best_paths = vector<bool>(num_vertices * collapsed_best_paths.size(), false);
    
    for (ushort path_idx = 0; path_idx < collapsed_best_paths.size(); path_idx++) {

        for (auto & vertex: collapsed_best_paths.at(path_idx)->getPath()) {

            assert(!(best_paths.at((path_idx * num_vertices) + vertex.first)));
            best_paths.at((path_idx * num_vertices) + vertex.first) = true;
        }

        delete collapsed_best_paths.at(path_idx);
    }
}

template <uchar kmer_size>
void VariantClusterKmerGraph<kmer_size>::mergePaths(vector<VariantClusterGraphPath<kmer_size> *> * main_paths, vector<VariantClusterGraphPath<kmer_size> *> * input_paths, const bool is_last_vertex_use) {

    uint main_path_size = main_paths->size();
    main_paths->reserve(main_path_size + input_paths->size());
 
    for (auto & input_path: *input_paths) {

        bool has_equal = false;

        for (uint main_path_idx = 0; main_path_idx < main_path_size; main_path_idx++) {

            if (*(main_paths->at(main_path_idx)) == *input_path) {

                if (swapRedundantPath(*(main_paths->at(main_path_idx)), *input_path)) {

                    delete main_paths->at(main_path_idx);

                    if (is_last_vertex_use) {

                        main_paths->at(main_path_idx) = input_path;

                    } else {

                        main_paths->at(main_path_idx) = new VariantClusterGraphPath<kmer_size>(*input_path);                        
                    }
                
                } else if (is_last_vertex_use) {

                    delete input_path;
                }
                
                has_redundant_sequence = true;
                has_equal = true;
                break;
            }
        }

        if (!has_equal) {

            if (is_last_vertex_use) {

                main_paths->emplace_back(input_path);

            } else {

                main_paths->emplace_back(new VariantClusterGraphPath<kmer_size>(*input_path));
            }
        }
    } 
}

template <uchar kmer_size>
bool VariantClusterKmerGraph<kmer_size>::swapRedundantPath(const VariantClusterGraphPath<kmer_size> & redundant_path_1, const VariantClusterGraphPath<kmer_size> & redundant_path_2) {

    if (redundant_path_1.getPath().size() > redundant_path_2.getPath().size()) {

        return false;
    
    } else if (redundant_path_1.getPath().size() < redundant_path_2.getPath().size()) {

        return true;
    
    } else {

        for (uint vertex_idx = 0; vertex_idx < redundant_path_1.getPath().size(); vertex_idx++) {

            if (redundant_path_1.getPath().at(vertex_idx).first < redundant_path_2.getPath().at(vertex_idx).first) {

                return false;
            
            } else if (redundant_path_1.getPath().at(vertex_idx).first > redundant_path_2.getPath().at(vertex_idx).first) {

                return true;
            }
        }
    }

    assert(false);
    return false;
}

template <uchar kmer_size>
void VariantClusterKmerGraph<kmer_size>::filterBestPaths(vector<VariantClusterGraphPath<kmer_size> *> * best_paths, const uint max_paths, const bool force_sort) {

    assert(!(best_paths->empty()));

    if ((best_paths->size() > max_paths) or force_sort) {

        bool is_first_pass = true;
        unordered_set<uint> covered_vertices;

        auto sorted_best_paths_end_it = best_paths->begin();

        while ((sorted_best_paths_end_it != best_paths->end()) and ((sorted_best_paths_end_it - best_paths->begin()) < max_paths)) {

            auto cur_best_path_it = sorted_best_paths_end_it;
            auto cur_best_path_kmer_score = (*cur_best_path_it)->getPathKmerScore();
            auto cur_best_path_vertex_score = (*cur_best_path_it)->getPathVertexScore(covered_vertices);
            
            auto best_paths_it = sorted_best_paths_end_it;
            best_paths_it++;

            while (best_paths_it != best_paths->end()) {

                auto cur_path_kmer_score = (*best_paths_it)->getPathKmerScore();
                auto cur_path_vertex_score = (*best_paths_it)->getPathVertexScore(covered_vertices);

                if ((cur_path_vertex_score > 0) or !is_first_pass) {

                    if ((Utils::doubleCompare(cur_path_kmer_score, cur_best_path_kmer_score) and (cur_path_vertex_score > cur_best_path_vertex_score)) or (cur_path_kmer_score > cur_best_path_kmer_score)) {

                        cur_best_path_it = best_paths_it;         
                        cur_best_path_kmer_score = cur_path_kmer_score;
                        cur_best_path_vertex_score = cur_path_vertex_score;
                    }
                }

                best_paths_it++;
            }

            if (cur_best_path_vertex_score > 0) {

                (*cur_best_path_it)->updateCoveredVertices(&covered_vertices);
            }

            if (*sorted_best_paths_end_it != *cur_best_path_it) {

                swap(*sorted_best_paths_end_it, *cur_best_path_it);
            }

            if (is_first_pass and (cur_best_path_vertex_score == 0)) {

                is_first_pass = false;
            
            } else {

                sorted_best_paths_end_it++;
            }
        }
        
        assert(covered_vertices.size() <= boost::num_vertices(graph));
        
        auto best_paths_it = best_paths->begin();
        advance(best_paths_it, min(static_cast<uint>(best_paths->size()), max_paths));

        while (best_paths_it != best_paths->end()) {
            
            delete *best_paths_it;
            best_paths_it++;
        }

        best_paths->resize(min(static_cast<uint>(best_paths->size()), max_paths));
    }
}

template <uchar kmer_size>
vector<VariantClusterGraphPath<kmer_size> *> VariantClusterKmerGraph<kmer_size>::collapsePaths(vector<vector<VariantClusterGraphPath<kmer_size> *> > * best_paths_per_sample, mt19937 * prng, const ushort max_sample_paths) {

    vector<VariantClusterGraphPath<kmer_size> *> collapsed_best_paths;
    collapsed_best_paths.reserve(max_total_paths);

    shuffle(best_paths_per_sample->begin(), best_paths_per_sample->end(), *prng);

    for (ushort sample_path_idx = 0; sample_path_idx < max_sample_paths; sample_path_idx++) {

        for (auto & sample_paths: *best_paths_per_sample) {

            assert(sample_paths.size() <= max_sample_paths);

            if (sample_path_idx >= sample_paths.size()) {

                continue;
            }

            assert(sample_paths.at(sample_path_idx));

            if (collapsed_best_paths.size() < max_total_paths) {

                bool has_equal = false;
                auto collapsed_best_paths_it = collapsed_best_paths.begin();

                while (collapsed_best_paths_it != collapsed_best_paths.end()) {

                    if (sample_paths.at(sample_path_idx)->getPath() == (*collapsed_best_paths_it)->getPath()) {

                        delete sample_paths.at(sample_path_idx);

                        has_equal = true; 
                        break;

                    } else if (*(sample_paths.at(sample_path_idx)) == **collapsed_best_paths_it) {

                        if (swapRedundantPath(**collapsed_best_paths_it, *(sample_paths.at(sample_path_idx)))) {

                            delete *collapsed_best_paths_it;
                            *collapsed_best_paths_it = sample_paths.at(sample_path_idx);
                        
                        } else {

                            delete sample_paths.at(sample_path_idx);
                        }

                        has_equal = true;
                        has_redundant_sequence = true;
                        break;
                    }

                    collapsed_best_paths_it++;
                }

                if (!has_equal) {

                    collapsed_best_paths.push_back(sample_paths.at(sample_path_idx));
                } 

            } else {

                delete sample_paths.at(sample_path_idx);
            }
        } 
    }

    return collapsed_best_paths; 
} 

template <uchar kmer_size>
typename VariantClusterKmerGraph<kmer_size>::KmerMultiplicitiesIndex VariantClusterKmerGraph<kmer_size>::indexKmerMultiplicities(KmerHash * kmer_hash) {

    KmerMultiplicitiesIndex kmer_multiplicities_index;
    
    KmerPair<kmer_size> kmer_pair;
    bitset<2> nt_bits;

    const uint num_vertices = boost::num_vertices(graph);

    assert(!(best_paths.empty()));
    assert((best_paths.size() % num_vertices) == 0);

    const ushort num_best_paths = best_paths.size() / num_vertices;

    for (ushort path_idx = 0; path_idx < num_best_paths; path_idx++) {

        kmer_pair.reset();

        for (uint vertex_idx = 0; vertex_idx < num_vertices; vertex_idx++) {

            if (best_paths.at((path_idx * num_vertices) + vertex_idx)) {

                if (graph[vertex_idx].is_disconnected) {

                    kmer_pair.reset();
                }

                assert((graph[vertex_idx].sequence.size() % 2) == 0);
                auto sequence_it = graph[vertex_idx].sequence.begin();

                while (sequence_it != graph[vertex_idx].sequence.end()) {

                    nt_bits.set(0, *sequence_it);
                    sequence_it++;

                    nt_bits.set(1, *sequence_it);
                    sequence_it++;

                    if (kmer_pair.move(nt_bits)) {     

                        auto kmer_count = static_cast<BasicKmerHash<kmer_size> *>(kmer_hash)->findKmer(kmer_pair.getLexicographicalLowestKmer());

                        if (kmer_count) {

                            auto kmer_multiplicities_index_it = kmer_multiplicities_index.index.insert(kmer_pair.getLexicographicalLowestKmer(), KmerPathInfo(kmer_count, num_best_paths), true);

                            if (kmer_multiplicities_index_it.second) {

                                kmer_multiplicities_index.num_kmers++;
                            }

                            if ((*kmer_multiplicities_index_it.first).second.multiplicities.at(path_idx) < Utils::uchar_overflow) {

                                (*kmer_multiplicities_index_it.first).second.multiplicities.at(path_idx)++;
                            } 
                        }
                    }
                }
            }
        }
    }

    return kmer_multiplicities_index;
}

template <uchar kmer_size>
typename VariantClusterKmerGraph<kmer_size>::VariantKmerMultiplicitiesIndex VariantClusterKmerGraph<kmer_size>::indexVariantKmerMultiplicities(KmerHash * kmer_hash) {

    VariantKmerMultiplicitiesIndex variant_kmer_multiplicities_index;
    
    KmerPair<kmer_size> kmer_pair;
    bitset<2> nt_bits;

    const uint num_vertices = boost::num_vertices(graph);

    assert(!(best_paths.empty()));
    assert((best_paths.size() % num_vertices) == 0);

    const ushort num_best_paths = best_paths.size() / num_vertices;

    unordered_map<ushort, pair<bool, uint> > running_variants;

    for (ushort path_idx = 0; path_idx < num_best_paths; path_idx++) {

        kmer_pair.reset();

        uint num_nucleotides = 0;
        running_variants.clear();

        for (uint vertex_idx = 0; vertex_idx < num_vertices; vertex_idx++) {

            if (best_paths.at((path_idx * num_vertices) + vertex_idx)) {

                if (graph[vertex_idx].variant_allele_idx.first != Utils::ushort_overflow) {

                    assert(graph[vertex_idx].variant_allele_idx.second != Utils::ushort_overflow);
                    
                    auto running_variants_it = running_variants.emplace(graph[vertex_idx].variant_allele_idx.first, make_pair(graph[vertex_idx].variant_allele_idx.second == 0, num_nucleotides + kmer_size - 1));
                    assert((running_variants_it.first->second.second - kmer_size + 1) == num_nucleotides);

                    running_variants_it.first->second.second += (graph[vertex_idx].sequence.size() / 2);
                }

                for (auto & reference_variant_idx: graph[vertex_idx].reference_variant_indices) {

                    assert(reference_variant_idx != Utils::ushort_overflow);

                    auto running_variants_it = running_variants.find(reference_variant_idx);

                    if (running_variants_it != running_variants.end()) {

                        if (running_variants_it->second.first) {

                            assert((running_variants_it->second.second - kmer_size + 1) == num_nucleotides);
                            running_variants_it->second.second += (graph[vertex_idx].sequence.size() / 2);
                        }
                    }
                }

                if (graph[vertex_idx].is_disconnected) {

                    kmer_pair.reset();
                }

                assert((graph[vertex_idx].sequence.size() % 2) == 0);
                auto sequence_it = graph[vertex_idx].sequence.begin();

                while (sequence_it != graph[vertex_idx].sequence.end()) {

                    nt_bits.set(0, *sequence_it);
                    sequence_it++;

                    nt_bits.set(1, *sequence_it);
                    sequence_it++;

                    if (kmer_pair.move(nt_bits)) {     

                        auto kmer_counts = static_cast<BasicKmerHash<kmer_size> *>(kmer_hash)->findKmer(kmer_pair.getLexicographicalLowestKmer());

                        bool add_kmer = true;

                        if (kmer_counts) {

                            assert(kmer_counts->hasClusterOccurrence());

                            if (kmer_counts->isExcluded()) {

                                has_excluded_kmer = true;
                                add_kmer = false;
                            } 
                        }

                        if (add_kmer) {

                            auto variant_kmer_multiplicities_index_it = variant_kmer_multiplicities_index.index.insert(kmer_pair.getLexicographicalLowestKmer(), VariantKmerPathInfo(kmer_counts, num_best_paths), true);

                            if (variant_kmer_multiplicities_index_it.second) {
                        
                                if (kmer_counts) {
                                
                                    if (kmer_counts->hasMulticlusterOccurrence()) {

                                        variant_kmer_multiplicities_index.num_multicluster_kmers++;
                                    
                                    } else {

                                        variant_kmer_multiplicities_index.num_unique_kmers++;
                                    }
                                
                                } else {

                                    variant_kmer_multiplicities_index.num_unique_kmers++;                        
                                }
                            }

                            if ((*variant_kmer_multiplicities_index_it.first).second.multiplicities.at(path_idx) < Utils::uchar_overflow) {

                                (*variant_kmer_multiplicities_index_it.first).second.multiplicities.at(path_idx)++;
                            }

                            auto running_variants_it = running_variants.begin();

                            while (running_variants_it != running_variants.end()) {

                                if (running_variants_it->second.second <= num_nucleotides) {

                                    running_variants_it = running_variants.erase(running_variants_it);
                                    continue;
                                }

                                if ((*variant_kmer_multiplicities_index_it.first).second.variant_path_indices.empty()) {

                                    (*variant_kmer_multiplicities_index_it.first).second.variant_path_indices.emplace_back(running_variants_it->first, vector<bool>(num_best_paths, false));
                                    (*variant_kmer_multiplicities_index_it.first).second.variant_path_indices.back().second.at(path_idx) = true;

                                } else {

                                    bool new_variant_idx = true;

                                    for (auto & variant_path_idx: (*variant_kmer_multiplicities_index_it.first).second.variant_path_indices) {

                                        if (variant_path_idx.first == running_variants_it->first) {

                                            assert(variant_path_idx.second.size() == num_best_paths);
                                            variant_path_idx.second.at(path_idx) = true;

                                            new_variant_idx = false;
                                            break;
                                        }
                                    }

                                    if (new_variant_idx) {

                                        (*variant_kmer_multiplicities_index_it.first).second.variant_path_indices.emplace_back(running_variants_it->first, vector<bool>(num_best_paths, false));
                                        (*variant_kmer_multiplicities_index_it.first).second.variant_path_indices.back().second.at(path_idx) = true;
                                    }
                                }

                                running_variants_it++;
                            }
                        }
                    }

                    num_nucleotides++;
                }
            }
        }
    }

    return variant_kmer_multiplicities_index;
}


bool VariantClusterGraphCompare(VariantClusterGraph * first_variant_cluster_graph, VariantClusterGraph * second_variant_cluster_graph) { 

    assert(first_variant_cluster_graph);
    assert(second_variant_cluster_graph);

    return (first_variant_cluster_graph->getInfo().size() > second_variant_cluster_graph->getInfo().size());
}

template class VariantClusterKmerGraph<31>;
template class VariantClusterKmerGraph<39>;
template class VariantClusterKmerGraph<47>;
template class VariantClusterKmerGraph<55>;
template class VariantClusterKmerGraph<63>;


