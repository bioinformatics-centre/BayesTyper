
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

static const uchar max_num_alternative_alleles = 14;
static const ushort haplotype_candidate_nbest_multiplicity = 3;


VariantClusterGraph::VariantClusterGraph(VariantCluster * variant_cluster, const string & chromosome_sequence, const uchar kmer_size) {

    has_ambiguous_nucleotide = false;
    has_complex_region = false;
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

    while (vit.first != vit.second) {

        KmerPair<Utils::small_kmer_size> cur_pair;
        num_unique_smallmers += countVertexSmallmers(smallmer_set, cur_pair, Utils::uint_overflow, *vit.first, 0);

        vit.first++;
    }

    return num_unique_smallmers;
}

ulong VariantClusterGraph::countVertexSmallmers(Utils::SmallmerSet * smallmer_set, KmerPair<Utils::small_kmer_size> cur_pair, uint remaining_prefix, const vertex_t cur_vertex, uchar num_alternative_alleles) {

    if ((graph[cur_vertex].variant_allele_idx.first != Utils::ushort_overflow) and (graph[cur_vertex].variant_allele_idx.second != 0)) {

        assert(graph[cur_vertex].variant_allele_idx.second != Utils::ushort_overflow);
        num_alternative_alleles++;
    }

    ulong num_unique_smallmers = 0;

    if (max_num_alternative_alleles < num_alternative_alleles) {

        has_complex_region = true;
        return num_unique_smallmers;
    }

    if (graph[cur_vertex].is_disconnected) {

        cur_pair.reset();
    }

    bitset<2> nt_bits;

    assert((graph[cur_vertex].sequence.size() % 2) == 0);    
    auto sequence_it = graph[cur_vertex].sequence.cbegin();

    while (sequence_it != graph[cur_vertex].sequence.cend()) {

        nt_bits.set(0, *sequence_it);
        sequence_it++;

        nt_bits.set(1, *sequence_it);
        sequence_it++;

       if (cur_pair.move(nt_bits)) { 

            num_unique_smallmers += smallmer_set->insert(cur_pair.getForwardKmer());
            num_unique_smallmers += smallmer_set->insert(cur_pair.getReverseComplementKmer());
        } 

        remaining_prefix--;

        if (remaining_prefix == 0) {

            return num_unique_smallmers;
        } 
    }

    auto out_eit = boost::out_edges(cur_vertex, graph);

    while (out_eit.first != out_eit.second) {
                    
        if (remaining_prefix < Utils::small_kmer_size) {

            num_unique_smallmers += countVertexSmallmers(smallmer_set, cur_pair, remaining_prefix, boost::target(*out_eit.first, graph), num_alternative_alleles);

        } else {

            num_unique_smallmers += countVertexSmallmers(smallmer_set, cur_pair, Utils::small_kmer_size - 1, boost::target(*out_eit.first, graph), num_alternative_alleles);
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
double VariantClusterKmerGraph<kmer_size>::countKmers(KmerHash * kmer_hash, const uint variant_cluster_group_idx, const uint prng_seed, const ushort num_samples, const ushort num_haplotype_candidates_per_sample) {
    
    mt19937 prng = mt19937(prng_seed);

    auto best_paths = findBestPaths(kmer_hash, &prng, num_samples, num_haplotype_candidates_per_sample);   

    assert(!(best_paths.empty()));
    assert(best_paths.size() <= static_cast<uint>(num_samples * num_haplotype_candidates_per_sample));

    auto kmer_multiplicities_index = indexKmerMultiplicities(kmer_hash, best_paths);

    for (auto & path: best_paths) {

        delete path;
    }

    for (auto &kit: kmer_multiplicities_index.index) {

        assert(kit.second.kmer_counts);
        assert(kit.second.multiplicities.size() == best_paths.size());

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

    return (kmer_multiplicities_index.num_kmers * (((best_paths.size() + 1) * best_paths.size()) / 2 + best_paths.size()));
}

template <uchar kmer_size>
void VariantClusterKmerGraph<kmer_size>::getBestHaplotypeCandidates(KmerHash * kmer_hash, VariantClusterHaplotypes * variant_cluster_haplotypes, const uint prng_seed, const ushort num_samples, const ushort num_haplotype_candidates_per_sample, const uchar num_genomic_rate_gc_bias_bins) {

    mt19937 prng = mt19937(prng_seed);

    auto best_paths = findBestPaths(kmer_hash, &prng, num_samples, num_haplotype_candidates_per_sample);   

    assert(!(best_paths.empty()));
    assert(best_paths.size() <= static_cast<uint>(num_samples * num_haplotype_candidates_per_sample));
    
    assert(variant_cluster_haplotypes->empty());

    variant_cluster_haplotypes->haplotypes.reserve(best_paths.size());

    for (ushort path_idx = 0; path_idx < best_paths.size(); path_idx++) {

        variant_cluster_haplotypes->haplotypes.emplace_back(variant_cluster_info.size());

        for (auto & vertex: best_paths.at(path_idx)->getPath()) {

            if (vertex->variant_allele_idx.first != Utils::ushort_overflow) {

                assert(vertex->variant_allele_idx.second != Utils::ushort_overflow);

                if (!(vertex->is_disconnected)) {

                    assert(variant_cluster_haplotypes->haplotypes.back().variant_allele_indices.at(vertex->variant_allele_idx.first) == Utils::ushort_overflow);
                    variant_cluster_haplotypes->haplotypes.back().variant_allele_indices.at(vertex->variant_allele_idx.first) = vertex->variant_allele_idx.second;
                }

                assert(variant_cluster_haplotypes->haplotypes.back().variant_allele_indices.at(vertex->variant_allele_idx.first) == vertex->variant_allele_idx.second);
            }

            if (vertex->nested_variant_cluster_index != Utils::uint_overflow) {

                assert(vertex->is_disconnected);
                variant_cluster_haplotypes->haplotypes.back().nested_variant_cluster_indices.push_back(vertex->nested_variant_cluster_index);
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

    auto variant_kmer_multiplicities_index = indexVariantKmerMultiplicities(kmer_hash, best_paths);

    for (auto & path: best_paths) {

        delete path;
    }

    variant_cluster_haplotypes->haplotype_unique_kmer_multiplicities = Eigen::MatrixXuchar(variant_kmer_multiplicities_index.num_unique_kmers, best_paths.size()); 
    variant_cluster_haplotypes->haplotype_multicluster_kmer_multiplicities = Eigen::MatrixXuchar(variant_kmer_multiplicities_index.num_multicluster_kmers, best_paths.size()); 

    variant_cluster_haplotypes->unique_kmers.reserve(variant_kmer_multiplicities_index.num_unique_kmers);
    variant_cluster_haplotypes->multicluster_kmers.reserve(variant_kmer_multiplicities_index.num_multicluster_kmers);

    uint unique_row_idx = 0;
    uint multicluster_row_idx = 0;

    for (auto &kit: variant_kmer_multiplicities_index.index) {

        assert(kit.second.multiplicities.size() == best_paths.size());
        
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
vector<VariantClusterGraphPath<kmer_size> *> VariantClusterKmerGraph<kmer_size>::findBestPaths(KmerHash * kmer_hash, mt19937 * prng, const ushort num_samples, const ushort num_haplotype_candidates_per_sample) {

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

                path->addVertex(graph[*vit.first]);
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

            filterBestPaths(&(best_paths_per_sample->at(sample_idx)), static_cast<uint>(num_haplotype_candidates_per_sample * haplotype_candidate_nbest_multiplicity));
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

        filterBestPaths(&(best_paths_per_sample->at(sample_idx)), num_haplotype_candidates_per_sample);
        assert(!(best_paths_per_sample->at(sample_idx).empty()));            
    }

    return collapsePaths(best_paths_per_sample);
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
bool VariantClusterKmerGraph<kmer_size>::swapRedundantPath(VariantClusterGraphPath<kmer_size> & redundant_path_1, VariantClusterGraphPath<kmer_size> & redundant_path_2) {

    if (redundant_path_1.getPath().size() > redundant_path_2.getPath().size()) {

        return false;
    
    } else if (redundant_path_1.getPath().size() < redundant_path_2.getPath().size()) {

        return true;
    
    } else {

        for (uint vertex_idx = 0; vertex_idx < redundant_path_1.getPath().size(); vertex_idx++) {

            ushort allele_idx_1 = Utils::ushort_overflow;
            ushort allele_idx_2 = Utils::ushort_overflow;

            if (redundant_path_1.getPath().at(vertex_idx)->variant_allele_idx.first != Utils::ushort_overflow) {

                assert(allele_idx_1 == Utils::ushort_overflow);                
                
                allele_idx_1 = redundant_path_1.getPath().at(vertex_idx)->variant_allele_idx.second;
                assert(allele_idx_1 != Utils::ushort_overflow);
            }

            if (redundant_path_2.getPath().at(vertex_idx)->variant_allele_idx.first != Utils::ushort_overflow) {

                assert(allele_idx_2 == Utils::ushort_overflow);                
                
                allele_idx_2 = redundant_path_2.getPath().at(vertex_idx)->variant_allele_idx.second;
                assert(allele_idx_2 != Utils::ushort_overflow);
            }

            if (allele_idx_1 < allele_idx_2) {

                return false;
            
            } else if (allele_idx_1 > allele_idx_2) {

                return true;
            }
        }
    }

    assert(false);
    return false;
}

template <uchar kmer_size>
void VariantClusterKmerGraph<kmer_size>::filterBestPaths(vector<VariantClusterGraphPath<kmer_size> *> * best_paths, const uint max_paths) {

    assert(!(best_paths->empty()));

    if (best_paths->size() > max_paths) {

        list<pair<VariantClusterGraphPath<kmer_size> *, bool> > remaining_paths;

        for (auto & path: *best_paths) {

            remaining_paths.emplace_back(path, true);
        }

        best_paths->clear();

        bool is_first_pass = true;
        unordered_set<const VariantClusterGraphVertex *> covered_vertices;
       
        while (best_paths->size() < max_paths) {

            double best_path_kmer_score = 0;
            uint best_path_vertex_score = 0;

            typename list<pair<VariantClusterGraphPath<kmer_size> *, bool> >::iterator best_scoring_path_it = remaining_paths.end();

            auto remaining_paths_it = remaining_paths.begin();
            assert(remaining_paths_it != remaining_paths.end());

            while (remaining_paths_it != remaining_paths.end()) {

                double cur_path_kmer_score = remaining_paths_it->first->getPathKmerScore();
                uint cur_path_vertex_score = 0;

                if (remaining_paths_it->second) {

                    cur_path_vertex_score = remaining_paths_it->first->getPathVertexScore(covered_vertices);
                }

                if (cur_path_vertex_score == 0) {

                    remaining_paths_it->second = false;
                }

                assert(cur_path_kmer_score >= 0);
                assert(cur_path_vertex_score >= 0);

                if ((cur_path_vertex_score > 0) or !is_first_pass) {

                    if ((best_scoring_path_it == remaining_paths.end()) or (Utils::doubleCompare(cur_path_kmer_score, best_path_kmer_score) and (cur_path_vertex_score > best_path_vertex_score)) or (cur_path_kmer_score > best_path_kmer_score)) {

                        best_path_kmer_score = cur_path_kmer_score;
                        best_path_vertex_score = cur_path_vertex_score;
                        best_scoring_path_it = remaining_paths_it;         
                    }
                }

                remaining_paths_it++;
            }

            if (best_scoring_path_it == remaining_paths.end()) {

                assert(is_first_pass);
                is_first_pass = false;

            } else {

                assert(is_first_pass == (best_path_vertex_score > 0));

                if (is_first_pass) {

                    assert(best_scoring_path_it->second);
                    best_scoring_path_it->first->updateCoveredVertices(&covered_vertices);     
                }

                best_paths->push_back(best_scoring_path_it->first);
                remaining_paths.erase(best_scoring_path_it);
            }
        }

        assert(best_paths->size() == max_paths);

        assert(!(remaining_paths.empty()));        
        assert(covered_vertices.size() <= boost::num_vertices(graph));

        for (auto & path: remaining_paths) {
            
            delete path.first;
        }     
    }
}

template <uchar kmer_size>
vector<VariantClusterGraphPath<kmer_size> *> VariantClusterKmerGraph<kmer_size>::collapsePaths(vector<vector<VariantClusterGraphPath<kmer_size> *> > * best_paths_per_sample) {

    vector<VariantClusterGraphPath<kmer_size> *> collapsed_best_paths;

    for (auto &sample_paths: *best_paths_per_sample) {

        auto sample_paths_it = sample_paths.begin();

        while (sample_paths_it != sample_paths.end()) {

            assert(*sample_paths_it);
            bool has_equal = false;

            auto collapsed_best_paths_it = collapsed_best_paths.begin();

            while (collapsed_best_paths_it != collapsed_best_paths.end()) {

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

    return collapsed_best_paths; 
} 

template <uchar kmer_size>
typename VariantClusterKmerGraph<kmer_size>::KmerMultiplicitiesIndex VariantClusterKmerGraph<kmer_size>::indexKmerMultiplicities(KmerHash * kmer_hash, const vector<VariantClusterGraphPath<kmer_size> *> & collapsed_best_paths) {

    KmerMultiplicitiesIndex kmer_multiplicities_index;

    bitset<2> nt_bits;

    for (ushort path_idx = 0; path_idx < collapsed_best_paths.size(); path_idx++) {

        assert(collapsed_best_paths.at(path_idx));
        collapsed_best_paths.at(path_idx)->kmer_pair.reset();

        for (auto & vertex: collapsed_best_paths.at(path_idx)->getPath()) {

            if (vertex->is_disconnected) {

                collapsed_best_paths.at(path_idx)->kmer_pair.reset();
            }

            assert((vertex->sequence.size() % 2) == 0);
            auto sequence_it = vertex->sequence.begin();

            while (sequence_it != vertex->sequence.end()) {

                nt_bits.set(0, *sequence_it);
                sequence_it++;

                nt_bits.set(1, *sequence_it);
                sequence_it++;

                if (collapsed_best_paths.at(path_idx)->kmer_pair.move(nt_bits)) {     

                    auto kmer_count = static_cast<BasicKmerHash<kmer_size> *>(kmer_hash)->findKmer(collapsed_best_paths.at(path_idx)->kmer_pair.getLexicographicalLowestKmer());

                    if (kmer_count) {

                        auto kmer_multiplicities_index_it = kmer_multiplicities_index.index.insert(collapsed_best_paths.at(path_idx)->kmer_pair.getLexicographicalLowestKmer(), KmerPathInfo(kmer_count, collapsed_best_paths.size()), true);

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

    return kmer_multiplicities_index;
}

template <uchar kmer_size>
typename VariantClusterKmerGraph<kmer_size>::VariantKmerMultiplicitiesIndex VariantClusterKmerGraph<kmer_size>::indexVariantKmerMultiplicities(KmerHash * kmer_hash, const vector<VariantClusterGraphPath<kmer_size> *> & collapsed_best_paths) {

    VariantKmerMultiplicitiesIndex variant_kmer_multiplicities_index;
    
    bitset<2> nt_bits;
    unordered_map<ushort, uint> running_variants;

    for (ushort path_idx = 0; path_idx < collapsed_best_paths.size(); path_idx++) {

        uint num_nucleotides = 0;
        running_variants.clear();

        assert(collapsed_best_paths.at(path_idx));
        collapsed_best_paths.at(path_idx)->kmer_pair.reset();

        for (auto & vertex: collapsed_best_paths.at(path_idx)->getPath()) {

            if (vertex->variant_allele_idx.first != Utils::ushort_overflow) {

                assert(vertex->variant_allele_idx.second != Utils::ushort_overflow);
                
                auto running_variants_it = running_variants.emplace(vertex->variant_allele_idx.first, num_nucleotides + kmer_size - 1);
                assert((running_variants_it.first->second - kmer_size + 1) == num_nucleotides);

                running_variants_it.first->second += (vertex->sequence.size() / 2);
            }

            for (auto & reference_variant_idx: vertex->reference_variant_indices) {

                assert(reference_variant_idx != Utils::ushort_overflow);

                auto running_variants_it = running_variants.emplace(reference_variant_idx, num_nucleotides + kmer_size - 1);
                assert((running_variants_it.first->second - kmer_size + 1) == num_nucleotides);

                running_variants_it.first->second += (vertex->sequence.size() / 2);
            }

            if (vertex->is_disconnected) {

                collapsed_best_paths.at(path_idx)->kmer_pair.reset();
            }

            assert((vertex->sequence.size() % 2) == 0);
            auto sequence_it = vertex->sequence.begin();

            while (sequence_it != vertex->sequence.end()) {

                nt_bits.set(0, *sequence_it);
                sequence_it++;

                nt_bits.set(1, *sequence_it);
                sequence_it++;

                if (collapsed_best_paths.at(path_idx)->kmer_pair.move(nt_bits)) {     

                    auto kmer_counts = static_cast<BasicKmerHash<kmer_size> *>(kmer_hash)->findKmer(collapsed_best_paths.at(path_idx)->kmer_pair.getLexicographicalLowestKmer());

                    bool add_kmer = true;

                    if (kmer_counts) {

                        assert(kmer_counts->hasClusterOccurrence());

                        if (kmer_counts->isExcluded()) {

                            has_excluded_kmer = true;
                            add_kmer = false;
                        } 
                    }

                    if (add_kmer) {

                        auto variant_kmer_multiplicities_index_it = variant_kmer_multiplicities_index.index.insert(collapsed_best_paths.at(path_idx)->kmer_pair.getLexicographicalLowestKmer(), VariantKmerPathInfo(kmer_counts, collapsed_best_paths.size()), true);

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

                        assert((*variant_kmer_multiplicities_index_it.first).second.variant_path_indices.empty() == variant_kmer_multiplicities_index_it.second);

                        auto running_variants_it = running_variants.begin();

                        while (running_variants_it != running_variants.end()) {

                            if (running_variants_it->second <= num_nucleotides) {

                                running_variants_it = running_variants.erase(running_variants_it);
                                continue;
                            }

                            if ((*variant_kmer_multiplicities_index_it.first).second.variant_path_indices.empty()) {

                                (*variant_kmer_multiplicities_index_it.first).second.variant_path_indices.emplace_back(running_variants_it->first, vector<bool>(collapsed_best_paths.size(), false));
                                (*variant_kmer_multiplicities_index_it.first).second.variant_path_indices.back().second.at(path_idx) = true;

                            } else {

                                bool new_variant_idx = true;

                                for (auto & variant_path_idx: (*variant_kmer_multiplicities_index_it.first).second.variant_path_indices) {

                                    if (variant_path_idx.first == running_variants_it->first) {

                                        assert(variant_path_idx.second.size() == collapsed_best_paths.size());
                                        variant_path_idx.second.at(path_idx) = true;

                                        new_variant_idx = false;
                                        break;
                                    }
                                }

                                if (new_variant_idx) {

                                    (*variant_kmer_multiplicities_index_it.first).second.variant_path_indices.emplace_back(running_variants_it->first, vector<bool>(collapsed_best_paths.size(), false));
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

    return variant_kmer_multiplicities_index;
}


template class VariantClusterKmerGraph<31>;
template class VariantClusterKmerGraph<39>;
template class VariantClusterKmerGraph<47>;
template class VariantClusterKmerGraph<55>;
template class VariantClusterKmerGraph<63>;


