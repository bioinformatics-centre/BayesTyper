
/*
VariantClusterGraph.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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
#include "boost/functional/hash.hpp"

#include "VariantClusterGraph.hpp"
#include "Utils.hpp"
#include "VariantCluster.hpp"
#include "KmerHash.hpp"
#include "KmerBloom.hpp"
#include "VariantClusterGraphPath.hpp"
#include "Kmer.hpp"
#include "Nucleotide.hpp"
#include "KmerCounts.hpp"
#include "VariantClusterHaplotypes.hpp"
#include "VariantInfo.hpp"
#include "VariantClusterGraphVertex.hpp"


static const ushort min_num_sample_paths = 1;

VariantClusterGraph::VariantClusterGraph(VariantCluster * variant_cluster, const string & chrom_sequence) {

    num_path_kmers = 0;
    has_excluded_kmers = false;

    assert(variant_cluster->variants.size() < Utils::ushort_overflow);
    
    assert(variant_cluster_info.empty());
    variant_cluster_info.reserve(variant_cluster->variants.size());

    map<uint, pair<vector<vertex_t>, vector<ushort> > > added_vertices;
    unordered_set<ushort> reference_variant_indices;

    auto variants_it = variant_cluster->variants.begin();
    auto chrom_sequence_it = chrom_sequence.begin();

    vertex_t cur_vertex = boost::add_vertex(graph);    

    addVertices(&cur_vertex, vector<StringItPair>(1, StringItPair(chrom_sequence_it + variants_it->first - (Utils::kmer_size - 1), chrom_sequence_it + variants_it->first)), make_pair(Utils::ushort_overflow, Utils::ushort_overflow), reference_variant_indices, vector<uint>(), false);

    vertex_t prev_vertex = cur_vertex;

    assert(added_vertices.insert({variants_it->first, make_pair(vector<vertex_t>(1, cur_vertex), vector<ushort>())}).second);

    uint cur_last_position = 0;
    uint next_position = 0;

    ushort variant_counter = 0;

    while (variants_it != variant_cluster->variants.end()) {

        variant_cluster_info.emplace_back(variants_it->first + 1, variants_it->second.id, variants_it->second.has_dependency, variants_it->second.alt_alleles);

        const bool is_first_nucleotides_redundant = variants_it->second.num_redundant_nucleotides > 0;
        uint max_reference_length = 0;

        for (ushort alt_allele_idx = 0; alt_allele_idx < variants_it->second.alt_alleles.size(); alt_allele_idx++) {

            max_reference_length = max(max_reference_length, variants_it->second.alt_alleles.at(alt_allele_idx).ref_length);

            vertex_t next_vertex = boost::add_vertex(graph);
            assert(boost::add_edge(cur_vertex, next_vertex, graph).second);

            assert(variants_it->second.num_redundant_nucleotides <= variants_it->second.alt_alleles.at(alt_allele_idx).ref_length);
            assert(variants_it->second.num_redundant_nucleotides <= variants_it->second.alt_alleles.at(alt_allele_idx).sequence.size());

            addVertices(&next_vertex, vector<StringItPair>(1, StringItPair(variants_it->second.alt_alleles.at(alt_allele_idx).sequence.begin(), variants_it->second.alt_alleles.at(alt_allele_idx).sequence.end())), make_pair(variant_counter, alt_allele_idx + 1), reference_variant_indices, vector<uint>(), is_first_nucleotides_redundant);

            auto added_vertices_insert = added_vertices.insert({variants_it->first + variants_it->second.alt_alleles.at(alt_allele_idx).ref_length, make_pair(vector<vertex_t>(), vector<ushort>())});
            added_vertices_insert.first->second.first.push_back(next_vertex);
        }

        assert(max_reference_length > 0);

        auto added_vertices_it = added_vertices.find(variants_it->first + max_reference_length);
        assert(added_vertices_it != added_vertices.end());

        added_vertices_it->second.second.push_back(variant_counter);
        reference_variant_indices.insert(variant_counter);

        variants_it++;

        added_vertices_it = added_vertices.begin();
        bool more_edges = true;
        bool last_variant = false;

        if (variants_it != variant_cluster->variants.end()) {

            next_position = variants_it->first;    

        } else {

            next_position = Utils::uint_overflow;
            last_variant = true;
        }

        while (more_edges) {

            auto cur_position = added_vertices_it->first;
            auto next_vertices = added_vertices_it->second.first;

            for (auto & variant_idx: added_vertices_it->second.second) {
                
                reference_variant_indices.erase(variant_idx);
            }

            added_vertices.erase(added_vertices_it);

            if (added_vertices.empty()) {

                more_edges = false;

                if (last_variant) {

                    cur_last_position = cur_position + Utils::kmer_size - 1;

                } else {

                    cur_last_position = next_position;
                }                   

            } else {

                added_vertices_it = added_vertices.begin();
                cur_last_position = added_vertices_it->first;

                if (!last_variant and (cur_last_position > next_position)) {

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
                assert(contained_cluster_it->right_flank <= (cur_last_position - Utils::kmer_size));

                if (prev_contained_edge < Utils::uint_overflow) {

                    nested_variant_cluster_indices.emplace_back(prev_contained_edge);
                }

                contained_vertices.emplace_back(chrom_sequence_it + cur_position, chrom_sequence_it + contained_cluster_it->left_flank);

                prev_contained_edge = contained_cluster_it->cluster_idx;
                cur_position = contained_cluster_it->right_flank + 1;
                variant_cluster->contained_clusters.erase(contained_cluster_it++);
            }
 
            if (prev_contained_edge < Utils::uint_overflow) {

                nested_variant_cluster_indices.emplace_back(prev_contained_edge);
            }

            contained_vertices.emplace_back(chrom_sequence_it + cur_position, chrom_sequence_it + cur_last_position);
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

                addVertices(&cur_vertex, contained_vertices, pair<ushort, ushort>(variant_counter, 0), reference_variant_indices, nested_variant_cluster_indices, is_first_nucleotides_redundant);    

            } else {

                addVertices(&cur_vertex, contained_vertices, pair<ushort, ushort>(Utils::ushort_overflow, Utils::ushort_overflow), reference_variant_indices, nested_variant_cluster_indices, false);
            }

            auto added_vertices_insert = added_vertices.insert({cur_last_position, make_pair(vector<vertex_t>(), vector<ushort>())});
            added_vertices_insert.first->second.first.push_back(cur_vertex);
        }

        variant_counter++;
        prev_vertex = cur_vertex; 
    }

    // if (variant_cluster->left_flank == 13294446 - 1) {

    //     auto vit = boost::vertices(graph);

    //     cout << "\n\n### Vertices: \n" << endl;

    //     while (vit.first != vit.second) {

    //         cout << *vit.first << ": (" << graph[*vit.first].variant_allele_idx.first << "," << graph[*vit.first].variant_allele_idx.second << "), (";
            
    //         for (auto & reference_variant_index: graph[*vit.first].reference_variant_indices) {

    //             cout << reference_variant_index << ",";
    //         }

    //         cout << "), " << graph[*vit.first].nested_variant_cluster_index << ", " << graph[*vit.first].is_disconnected << ", " << graph[*vit.first].is_first_nucleotides_redundant << ", ";

    //         for (uint i = 0; i < graph[*vit.first].sequence.size(); i++) {

    //             cout << graph[*vit.first].sequence.at(i);
    //         }

    //         cout << endl;
    //         vit.first++;
    //     }

    //     cout << "\n### Edges: \n" << endl;

    //     auto eit = boost::edges(graph);

    //     while (eit.first != eit.second) {

    //         cout << boost::source(*eit.first, graph) << ">" << boost::target(*eit.first, graph) << std::endl;
    //         eit.first++;
    //     }

    //     cout << "\n" << endl;
    // }

    assert(variant_cluster->variants.size() == variant_counter);

    assert(added_vertices.size() == 1);
    assert(added_vertices.begin()->second.first.size() == 1);
    assert(added_vertices.begin()->second.second.empty());
}

void VariantClusterGraph::addVertices(vertex_t * cur_vertex, const vector<StringItPair> & vertex_sequences, const pair<ushort, ushort> & variant_allele_idx, const unordered_set<ushort> & reference_variant_indices, const vector<uint> & nested_variant_cluster_indices, const bool is_first_nucleotides_redundant) {

    assert(!vertex_sequences.empty());
    assert(vertex_sequences.size() == (nested_variant_cluster_indices.size() + 1));

    vector<ushort> vertex_reference_variant_indices;
    vertex_reference_variant_indices.reserve(reference_variant_indices.size());

    for (auto & reference_variant_idx: reference_variant_indices) {

        assert(reference_variant_idx != Utils::ushort_overflow);
    
        if (reference_variant_idx != variant_allele_idx.first) {

            vertex_reference_variant_indices.push_back(reference_variant_idx);
        }
    }

    initVertex(cur_vertex, vertex_sequences.front(), variant_allele_idx, vertex_reference_variant_indices, Utils::uint_overflow, is_first_nucleotides_redundant);
    
    for (uint vertex_idx = 1; vertex_idx < vertex_sequences.size(); vertex_idx++) {

        vertex_t prev_vertex = *cur_vertex;
        *cur_vertex = boost::add_vertex(graph);

        assert(boost::add_edge(prev_vertex, *cur_vertex, graph).second);

        initVertex(cur_vertex, vertex_sequences.at(vertex_idx), variant_allele_idx, vertex_reference_variant_indices, nested_variant_cluster_indices.at(vertex_idx - 1), false);
    }
}

void VariantClusterGraph::initVertex(vertex_t * cur_vertex, StringItPair vertex_sequence, const pair<ushort, ushort> & variant_allele_idx, const vector<ushort> & vertex_reference_variant_indices, const uint nested_variant_cluster_index, const bool is_first_nucleotides_redundant) {

    assert(static_cast<uint>(is_first_nucleotides_redundant) <= (vertex_sequence.second - vertex_sequence.first));

    graph[*cur_vertex].variant_allele_idx = variant_allele_idx;
    graph[*cur_vertex].reference_variant_indices = vertex_reference_variant_indices;
    graph[*cur_vertex].nested_variant_cluster_index = nested_variant_cluster_index;
    graph[*cur_vertex].is_first_nucleotides_redundant = is_first_nucleotides_redundant;

    if (nested_variant_cluster_index != Utils::uint_overflow) {

        assert((variant_allele_idx.second == 0) or (variant_allele_idx.second == Utils::ushort_overflow));
        graph[*cur_vertex].is_disconnected = true;

    } else {

        graph[*cur_vertex].is_disconnected = false;
    }

    graph[*cur_vertex].sequence.reserve((vertex_sequence.second - vertex_sequence.first) * 2);

    bool prev_is_disconnected = false;

    while (vertex_sequence.first != vertex_sequence.second) {

        auto nt_bits = Nucleotide::ntToBit<1>(*vertex_sequence.first);

        if (!nt_bits.second) {

            if (!prev_is_disconnected) {

                graph[*cur_vertex].sequence.shrink_to_fit();
                assert((graph[*cur_vertex].sequence.size() % 2) == 0);

                vertex_t prev_vertex = *cur_vertex;
                *cur_vertex = boost::add_vertex(graph);

                assert(boost::add_edge(prev_vertex, *cur_vertex, graph).second);

                graph[*cur_vertex].variant_allele_idx = variant_allele_idx;
                graph[*cur_vertex].reference_variant_indices = vertex_reference_variant_indices;
                graph[*cur_vertex].nested_variant_cluster_index = Utils::uint_overflow;
                graph[*cur_vertex].is_first_nucleotides_redundant = false;
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

    graph[*cur_vertex].sequence.shrink_to_fit();
    assert((graph[*cur_vertex].sequence.size() % 2) == 0);
}

const vector<VariantInfo> & VariantClusterGraph::getInfo() const {

    return variant_cluster_info;
}

bool VariantClusterGraph::hasExcludedKmers() {

    return has_excluded_kmers;
}

void VariantClusterGraph::findSamplePaths(KmerBloom<Utils::kmer_size> * sample_kmer_bloom, const uint prng_seed, const ushort max_sample_haplotypes) {
    
    mt19937 prng = mt19937(prng_seed);

    auto vit = boost::vertices(graph);

    auto eit_in = boost::in_edges(*vit.first, graph);
    auto eit_out = boost::out_edges(*vit.first, graph);

    assert(eit_in.first == eit_in.second);
    assert(eit_out.first != eit_out.second);

    unordered_map<vertex_t, vector<VariantClusterGraphPath *> > vertex_best_sample_paths;
    typename unordered_map<vertex_t, vector<VariantClusterGraphPath *> >::iterator vertex_best_sample_paths_it = vertex_best_sample_paths.end();

    unordered_map<vertex_t, vertex_t> visited_vertices;

    while (vit.first != vit.second) {

        auto vertex_best_sample_paths_emplace = vertex_best_sample_paths.emplace(*vit.first, vector<VariantClusterGraphPath *>());
        assert(vertex_best_sample_paths_emplace.second);

        vertex_best_sample_paths_it = vertex_best_sample_paths_emplace.first;

        eit_in = boost::in_edges(*vit.first, graph);

        if (eit_in.first == eit_in.second) {

            assert(*vit.first == 0);

            assert(graph[*vit.first].variant_allele_idx.first == Utils::ushort_overflow);
            assert(graph[*vit.first].variant_allele_idx.second == Utils::ushort_overflow);

            assert(graph[*vit.first].reference_variant_indices.empty());
            assert(graph[*vit.first].nested_variant_cluster_index == Utils::uint_overflow);

            vertex_best_sample_paths_it->second.emplace_back(new VariantClusterGraphPath(variant_cluster_info.size()));

        } else {
         
            while (eit_in.first != eit_in.second) {

                auto cur_source_vertex = boost::source(*eit_in.first, graph);

                assert(cur_source_vertex < *vit.first);

                auto visited_vertices_it = visited_vertices.find(cur_source_vertex);
            
                assert(visited_vertices_it != visited_vertices.end());
                assert(visited_vertices_it->first != visited_vertices_it->second);

                mergePaths(&(vertex_best_sample_paths_it->second), &(vertex_best_sample_paths.at(cur_source_vertex)), *vit.first == visited_vertices_it->second);

                eit_in.first++;
            }
        }

        vertex_best_sample_paths_it->second.shrink_to_fit();
        shuffle(vertex_best_sample_paths_it->second.begin(), vertex_best_sample_paths_it->second.end(), prng);

        for (auto & path: vertex_best_sample_paths_it->second) {

            path->addVertex(*vit.first, graph[*vit.first], sample_kmer_bloom);
        }

        filterPaths(&(vertex_best_sample_paths_it->second), max_sample_haplotypes, false);
        assert(!vertex_best_sample_paths_it->second.empty() and (vertex_best_sample_paths_it->second.size() <= max_sample_haplotypes));  

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

    assert(vertex_best_sample_paths_it != vertex_best_sample_paths.end());
    assert(vertex_best_sample_paths.size() == visited_vertices.size());
    assert(vertex_best_sample_paths.size() == boost::num_vertices(graph));

    filterPaths(&(vertex_best_sample_paths_it->second), max_sample_haplotypes, true);
    assert(!vertex_best_sample_paths_it->second.empty() and (vertex_best_sample_paths_it->second.size() <= max_sample_haplotypes));  

    addPathIndices(&(vertex_best_sample_paths_it->second));
}

void VariantClusterGraph::mergePaths(vector<VariantClusterGraphPath *> * main_paths, vector<VariantClusterGraphPath *> * input_paths, const bool is_last_edge) {

    const uint main_path_size = main_paths->size();
    main_paths->reserve(main_path_size + input_paths->size());
 
    for (auto & input_path: *input_paths) {

        bool is_redundant = false;

        for (uint main_path_idx = 0; main_path_idx < main_path_size; main_path_idx++) {

            if (isPathsRedundant(main_paths->at(main_path_idx)->getPath(), input_path->getPath())) {

                if (main_paths->at(main_path_idx)->getPath().size() < input_path->getPath().size()) {

                    delete main_paths->at(main_path_idx);

                    if (is_last_edge) {

                        main_paths->at(main_path_idx) = input_path;

                    } else {

                        main_paths->at(main_path_idx) = new VariantClusterGraphPath(*input_path);
                    }
                
                } else {

                    delete input_path;
                }

                is_redundant = true;
                break;
            }
        }

        if (!is_redundant) {

            if (is_last_edge) {

                main_paths->emplace_back(input_path);

            } else {

                main_paths->emplace_back(new VariantClusterGraphPath(*input_path));
            }
        } 
    } 
}

bool VariantClusterGraph::isPathsRedundant(const vector<typename VariantClusterGraphPath::PathVertexInfo> & path_1, const vector<typename VariantClusterGraphPath::PathVertexInfo> & path_2) const {

    auto path_it_1 = path_1.crbegin();
    auto path_it_2 = path_2.crbegin();

    assert(path_it_1 != path_1.crend());
    assert(path_it_2 != path_2.crend());

    auto sequence_it_1 = path_it_1->vertex->sequence.crbegin();
    auto sequence_it_2 = path_it_2->vertex->sequence.crbegin();

    bool is_disconnected_1 = false;
    bool is_disconnected_2 = false;

    while (true) {

        while (sequence_it_1 == path_it_1->vertex->sequence.crend()) {
            
            assert(path_it_1 != path_1.crend());

            if (path_it_1->vertex->is_disconnected) {

                is_disconnected_1 = true;
            }

            path_it_1++;

            if (path_it_1 != path_1.crend()) {

                sequence_it_1 = path_it_1->vertex->sequence.crbegin();
            
            } else {

                break;
            }
        }

        while (sequence_it_2 == path_it_2->vertex->sequence.crend()) {

            assert(path_it_2 != path_2.crend());

            if (path_it_2->vertex->is_disconnected) {

                is_disconnected_2 = true;
            }

            path_it_2++;

            if (path_it_2 != path_2.crend()) {

                sequence_it_2 = path_it_2->vertex->sequence.crbegin();
            
            } else {

                break;
            }
        }
        
        if (is_disconnected_1 != is_disconnected_2) {

            return false;
        }   

        is_disconnected_1 = false;
        is_disconnected_2 = false;

        if ((path_it_1 == path_1.crend()) or (path_it_2 == path_2.crend())) {

            break;
        }   

        if (sequence_it_1 == sequence_it_2) {

            sequence_it_1 = path_it_1->vertex->sequence.crend();
            sequence_it_2 = path_it_2->vertex->sequence.crend();
        }   

        while ((sequence_it_1 != path_it_1->vertex->sequence.crend()) and (sequence_it_2 != path_it_2->vertex->sequence.crend())) {

            if (*sequence_it_1 != *sequence_it_2) {

                return false;
            }

            sequence_it_1++;
            sequence_it_2++;
        }
    }

    if ((path_it_1 != path_1.crend()) or (path_it_2 != path_2.crend())) {

        return false;
    }

    return true;
}

void VariantClusterGraph::filterPaths(vector<VariantClusterGraphPath *> * paths, const uint max_paths, const bool is_complete) {

    assert(!paths->empty());

    if ((paths->size() > max_paths) or (is_complete and (paths->size() > min_num_sample_paths))) {

        bool is_first_pass = true;
        vector<bool> observed_covered_vertices(boost::num_vertices(graph), false);

        auto sorted_paths_end_it = paths->begin();

        while (sorted_paths_end_it != paths->end()) {

            auto cur_best_path_it = sorted_paths_end_it;
            auto cur_best_path_kmer_score = (*cur_best_path_it)->getKmerScore();
            auto cur_best_path_vertex_score = (*cur_best_path_it)->getVertexScore(observed_covered_vertices, is_complete);
            
            auto paths_it = sorted_paths_end_it;
            paths_it++;

            while (paths_it != paths->end()) {

                auto cur_path_kmer_score = (*paths_it)->getKmerScore();
                auto cur_path_vertex_score = (*paths_it)->getVertexScore(observed_covered_vertices, is_complete);

                if (is_first_pass) {

                    if (cur_path_vertex_score > 0) {

                        if ((Utils::doubleCompare(cur_path_kmer_score, cur_best_path_kmer_score) and (cur_path_vertex_score > cur_best_path_vertex_score)) or (cur_path_kmer_score > cur_best_path_kmer_score) or (cur_best_path_vertex_score == 0)) {

                            cur_best_path_it = paths_it;
                            cur_best_path_kmer_score = cur_path_kmer_score;
                            cur_best_path_vertex_score = cur_path_vertex_score;
                        } 
                    }

                } else if (!is_complete or (cur_path_vertex_score == (*paths_it)->getPath().size())) {

                    if (cur_path_kmer_score > cur_best_path_kmer_score) {

                        cur_best_path_it = paths_it;
                        cur_best_path_kmer_score = cur_path_kmer_score;
                        cur_best_path_vertex_score = cur_path_vertex_score;
                    }
                }

                paths_it++;
            }

            if (is_first_pass) {

                (*cur_best_path_it)->updateObservedCoveredVertices(&observed_covered_vertices, is_complete);

            } else if (is_complete and ((sorted_paths_end_it - paths->begin()) >= min_num_sample_paths) and (cur_best_path_vertex_score < (*cur_best_path_it)->getPath().size())) {

                break;
            }

            if (*sorted_paths_end_it != *cur_best_path_it) {

                swap(*sorted_paths_end_it, *cur_best_path_it);
            }

            if (is_first_pass and (cur_best_path_vertex_score == 0)) {

                is_first_pass = false;
                fill(observed_covered_vertices.begin(), observed_covered_vertices.end(), false);
            
            } else {

                sorted_paths_end_it++;

                if ((sorted_paths_end_it - paths->begin()) == max_paths) {

                    break;
                } 
            } 
        }

        const uint num_sorted_paths = sorted_paths_end_it - paths->begin();

        assert(num_sorted_paths >= min_num_sample_paths);
        assert(num_sorted_paths <= max_paths);
                
        while (sorted_paths_end_it != paths->end()) {
            
            delete *sorted_paths_end_it;
            sorted_paths_end_it++;
        }

        paths->resize(num_sorted_paths);
    }
}

void VariantClusterGraph::addPathIndices(vector<VariantClusterGraphPath *> * paths) {

    const uint num_vertices = boost::num_vertices(graph);

    assert((best_paths_indices.size() % num_vertices) == 0);

    const ushort num_best_paths = best_paths_indices.size() / num_vertices;
    
    best_paths_indices.reserve(best_paths_indices.size() + (num_vertices * paths->size()));

    vector<bool> redundant_paths(paths->size(), false);

    for (uint best_paths_idx = 0; best_paths_idx < num_best_paths; best_paths_idx++) {

        vector<typename VariantClusterGraphPath::PathVertexInfo> cur_best_path;
        cur_best_path.reserve(num_vertices / 2);

        for (uint vertex_idx = 0; vertex_idx < num_vertices; vertex_idx++) {

            if (best_paths_indices.at((best_paths_idx * num_vertices) + vertex_idx)) {

                cur_best_path.emplace_back(vertex_idx, &(graph[vertex_idx]));
            }
        }        
    
        for (ushort path_idx = 0; path_idx < paths->size(); path_idx++) {

            if (redundant_paths.at(path_idx)) {

                continue;
            }

            if (isPathsRedundant(cur_best_path, paths->at(path_idx)->getPath())) {

                if (cur_best_path.size() < paths->at(path_idx)->getPath().size()) {

                    vector<bool> path_indices(num_vertices, false);

                    for (auto & vertex: paths->at(path_idx)->getPath()) {

                        assert(!path_indices.at(vertex.index));
                        path_indices.at(vertex.index) = true;
                    }

                    copy(path_indices.begin(), path_indices.end(), best_paths_indices.begin() + best_paths_idx * num_vertices);
                }

                redundant_paths.at(path_idx) = true;
                break;
            }
        }
    }

    for (ushort path_idx = 0; path_idx < paths->size(); path_idx++) {

        if (!redundant_paths.at(path_idx)) {

            vector<bool> path_indices(num_vertices, false);

            for (auto & vertex: paths->at(path_idx)->getPath()) {

                assert(!path_indices.at(vertex.index));
                path_indices.at(vertex.index) = true;
            }

            best_paths_indices.insert(best_paths_indices.end(), path_indices.begin(), path_indices.end());
        }

        delete paths->at(path_idx);
    }

    best_paths_indices.shrink_to_fit();
} 

void VariantClusterGraph::countPathKmers(unordered_set<bitset<Utils::kmer_size * 2> > * path_kmers) {
    
    const uint num_vertices = boost::num_vertices(graph);

    assert(!best_paths_indices.empty());
    assert((best_paths_indices.size() % num_vertices) == 0);

    const ushort num_best_paths = best_paths_indices.size() / num_vertices;

    unordered_map<bitset<Utils::kmer_size * 2>, vector<uchar> > kmer_multiplicities;
    
    KmerPair<Utils::kmer_size> kmer_pair;
    bitset<2> nt_bits;

    for (ushort best_paths_idx = 0; best_paths_idx < num_best_paths; best_paths_idx++) {

        kmer_pair.reset();

        for (uint vertex_idx = 0; vertex_idx < num_vertices; vertex_idx++) {

            if (best_paths_indices.at((best_paths_idx * num_vertices) + vertex_idx)) {

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

                    if (kmer_pair.move(make_pair(nt_bits, true))) {

                        path_kmers->emplace(kmer_pair.getLexicographicalLowestKmer());
                    }
                }
            }
        }
    }
}

void VariantClusterGraph::classifyPathKmers(KmerCountsHash * kmer_hash, KmerBloom<Utils::kmer_size> * multigroup_kmer_bloom) {
    
    const uint num_vertices = boost::num_vertices(graph);

    assert(!best_paths_indices.empty());
    assert((best_paths_indices.size() % num_vertices) == 0);

    const ushort num_best_paths = best_paths_indices.size() / num_vertices;

    unordered_map<bitset<Utils::kmer_size * 2>, vector<uchar> > kmer_multiplicities;
    
    KmerPair<Utils::kmer_size> kmer_pair;
    bitset<2> nt_bits;

    for (ushort best_paths_idx = 0; best_paths_idx < num_best_paths; best_paths_idx++) {

        kmer_pair.reset();

        for (uint vertex_idx = 0; vertex_idx < num_vertices; vertex_idx++) {

            if (best_paths_indices.at((best_paths_idx * num_vertices) + vertex_idx)) {

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

                    if (kmer_pair.move(make_pair(nt_bits, true))) {
                    
                        auto kmer_multiplicities_it = kmer_multiplicities.emplace(kmer_pair.getLexicographicalLowestKmer(), vector<uchar>(num_best_paths, 0));

                        if (kmer_multiplicities_it.first->second.at(best_paths_idx) < Utils::uchar_overflow) {

                            kmer_multiplicities_it.first->second.at(best_paths_idx)++;
                        } 
                    }
                }
            }
        }
    }

    assert(num_path_kmers == 0);
    num_path_kmers += kmer_multiplicities.size();

    assert(has_excluded_kmers == false);

    auto kmer_multiplicities_it = kmer_multiplicities.begin();

    while (kmer_multiplicities_it != kmer_multiplicities.end()) {

        assert(kmer_multiplicities_it->second.size() == num_best_paths);

        uchar max_multiplicity = *max_element(kmer_multiplicities_it->second.begin(), kmer_multiplicities_it->second.end());
        assert(max_multiplicity > 0); 

        auto hash_lock = kmer_hash->getKmerLock(kmer_multiplicities_it->first);
        auto kmer_counts = kmer_hash->findKmer(kmer_multiplicities_it->first);    

        if (!kmer_counts and (max_multiplicity > Utils::bit7_overflow)) {

            auto kmer_added = kmer_hash->addKmer(kmer_multiplicities_it->first, true);
           
            assert(kmer_added.first);
            assert(kmer_added.second);

            kmer_counts = kmer_added.first;
        }

        if (kmer_counts) {

            kmer_counts->addClusterMultiplicity(max_multiplicity, multigroup_kmer_bloom->lookup(kmer_multiplicities_it->first));

            if (kmer_counts->isExcluded()) {

                has_excluded_kmers = true;
            }
        }

        kmer_multiplicities_it++;
    }
}

VariantClusterHaplotypes VariantClusterGraph::getHaplotypeCandidates(KmerCountsHash * kmer_hash, const uchar num_genomic_rate_gc_bias_bins) {

    assert(num_genomic_rate_gc_bias_bins == 1);

    VariantClusterHaplotypes variant_cluster_haplotypes;

    const uint num_vertices = boost::num_vertices(graph);

    assert(!best_paths_indices.empty());
    assert((best_paths_indices.size() % num_vertices) == 0);
    
    const ushort num_best_paths = best_paths_indices.size() / num_vertices;

    variant_cluster_haplotypes.haplotypes.reserve(best_paths_indices.size());

    variant_cluster_haplotypes.haplotype_kmer_multiplicities = Utils::MatrixXuchar::Zero(num_path_kmers, num_best_paths); 
    
    variant_cluster_haplotypes.kmers.reserve(num_path_kmers);
    variant_cluster_haplotypes.unique_kmer_indices.reserve(num_path_kmers);
    variant_cluster_haplotypes.multicluster_kmer_indices.reserve(num_path_kmers);

    unordered_map<bitset<Utils::kmer_size * 2>, uint> kmer_row_indices;

    KmerPair<Utils::kmer_size> kmer_pair;
    bitset<2> nt_bits;

    uint num_nucleotides = 0;
    unordered_map<pair<ushort, ushort>, pair<uint, uint>, boost::hash<pair<ushort, ushort> > > running_variants;

    for (ushort best_paths_idx = 0; best_paths_idx < num_best_paths; best_paths_idx++) {

        variant_cluster_haplotypes.haplotypes.emplace_back(variant_cluster_info.size());

        kmer_pair.reset();

        num_nucleotides = 0;
        running_variants.clear();

        for (uint vertex_idx = 0; vertex_idx < num_vertices; vertex_idx++) {

            if (best_paths_indices.at((best_paths_idx * num_vertices) + vertex_idx)) {

                if (graph[vertex_idx].variant_allele_idx.first != Utils::ushort_overflow) {

                    assert(graph[vertex_idx].variant_allele_idx.second != Utils::ushort_overflow);

                    if (!graph[vertex_idx].is_disconnected) {

                        assert(variant_cluster_haplotypes.haplotypes.back().variant_allele_indices.at(graph[vertex_idx].variant_allele_idx.first) == Utils::ushort_overflow);
                        variant_cluster_haplotypes.haplotypes.back().variant_allele_indices.at(graph[vertex_idx].variant_allele_idx.first) = graph[vertex_idx].variant_allele_idx.second;
                    }

                    assert(variant_cluster_haplotypes.haplotypes.back().variant_allele_indices.at(graph[vertex_idx].variant_allele_idx.first) == graph[vertex_idx].variant_allele_idx.second);
                    
                    auto running_variants_it = running_variants.emplace(graph[vertex_idx].variant_allele_idx, make_pair(num_nucleotides + static_cast<uint>(graph[vertex_idx].is_first_nucleotides_redundant), num_nucleotides + Utils::kmer_size - 1));
                    assert((running_variants_it.first->second.second - Utils::kmer_size + 1) == num_nucleotides);

                    running_variants_it.first->second.second += graph[vertex_idx].sequence.size() / 2;
                }

                for (auto & reference_variant_idx: graph[vertex_idx].reference_variant_indices) {

                    assert(reference_variant_idx != Utils::ushort_overflow);

                    auto running_variants_it = running_variants.find(make_pair(reference_variant_idx, 0));

                    if (running_variants_it != running_variants.end()) {

                        assert((running_variants_it->second.second - Utils::kmer_size + 1) == num_nucleotides);
                        running_variants_it->second.second += graph[vertex_idx].sequence.size() / 2;
                    }
                }

                if (graph[vertex_idx].nested_variant_cluster_index != Utils::uint_overflow) {

                    assert(graph[vertex_idx].is_disconnected);
                    variant_cluster_haplotypes.haplotypes.back().nested_variant_cluster_indices.push_back(graph[vertex_idx].nested_variant_cluster_index);
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

                    if (kmer_pair.move(make_pair(nt_bits, true))) {     

                        auto lowest_kmer = kmer_pair.getLexicographicalLowestKmer();
                        auto kmer_counts = kmer_hash->findKmer(lowest_kmer);

                        bool is_multicluster = false;

                        if (kmer_counts) {

                            assert(kmer_counts->hasClusterOccurrence());

                            if (kmer_counts->isExcluded()) {

                                num_nucleotides++;
                                continue;
                            }

                            is_multicluster = kmer_counts->hasMulticlusterOccurrence();
                        }

                        auto kmer_row_indices_it = kmer_row_indices.emplace(lowest_kmer, kmer_row_indices.size());
                        assert(kmer_row_indices_it.first->second < num_path_kmers);

                        variant_cluster_haplotypes.haplotype_kmer_multiplicities(kmer_row_indices_it.first->second, best_paths_idx)++;
                        assert(variant_cluster_haplotypes.haplotype_kmer_multiplicities(kmer_row_indices_it.first->second, best_paths_idx) <= Utils::bit7_overflow);

                        if (kmer_row_indices_it.second) {

                            variant_cluster_haplotypes.kmers.emplace_back(kmer_counts, Nucleotide::gcBiasBin<Utils::kmer_size>(lowest_kmer, num_genomic_rate_gc_bias_bins));

                            if (is_multicluster) {

                                variant_cluster_haplotypes.multicluster_kmer_indices.emplace_back(kmer_row_indices_it.first->second);

                            } else {

                                variant_cluster_haplotypes.unique_kmer_indices.emplace_back(kmer_row_indices_it.first->second);
                            }
                        }

                        updateVariantPathIndices(&(variant_cluster_haplotypes.kmers.at(kmer_row_indices_it.first->second).variant_haplotype_indices), &running_variants, num_nucleotides, num_best_paths, best_paths_idx);
                    }

                    num_nucleotides++;
                }
            }
        }

        sort(variant_cluster_haplotypes.haplotypes.back().nested_variant_cluster_indices.begin(), variant_cluster_haplotypes.haplotypes.back().nested_variant_cluster_indices.end());

        for (ushort variant_idx = 0; variant_idx < variant_cluster_info.size(); variant_idx++) {

            if (variant_cluster_haplotypes.haplotypes.back().variant_allele_indices.at(variant_idx) == Utils::ushort_overflow) {

                assert(variant_cluster_info.at(variant_idx).has_dependency);
                variant_cluster_haplotypes.haplotypes.back().variant_allele_indices.at(variant_idx) = variant_cluster_info.at(variant_idx).numberOfAlleles() - 1;
            }
        }  
    }

    if (kmer_row_indices.size() < static_cast<uint>(variant_cluster_haplotypes.haplotype_kmer_multiplicities.rows())) {

        Eigen::NoChange_t no_change_t;
        variant_cluster_haplotypes.haplotype_kmer_multiplicities.conservativeResize(kmer_row_indices.size(), no_change_t);
    }

    variant_cluster_haplotypes.kmers.shrink_to_fit();
    variant_cluster_haplotypes.unique_kmer_indices.shrink_to_fit();
    variant_cluster_haplotypes.multicluster_kmer_indices.shrink_to_fit();

    assert(kmer_row_indices.size() == static_cast<uint>(variant_cluster_haplotypes.haplotype_kmer_multiplicities.rows()));
    assert(kmer_row_indices.size() == variant_cluster_haplotypes.kmers.size());
    assert(kmer_row_indices.size() == (variant_cluster_haplotypes.unique_kmer_indices.size() + variant_cluster_haplotypes.multicluster_kmer_indices.size()));

    assert(variant_cluster_haplotypes.nested_variant_cluster_dependency.empty());

    for (uint vertex_idx = 0; vertex_idx < num_vertices; vertex_idx++) {

        if (graph[vertex_idx].nested_variant_cluster_index != Utils::uint_overflow) {

            auto nested_variant_cluster_dependency_it = variant_cluster_haplotypes.nested_variant_cluster_dependency.emplace(graph[vertex_idx].nested_variant_cluster_index, vector<ushort>());
            assert(nested_variant_cluster_dependency_it.second);

            if (graph[vertex_idx].variant_allele_idx.first != Utils::ushort_overflow) {

                nested_variant_cluster_dependency_it.first->second.push_back(graph[vertex_idx].variant_allele_idx.first);
            }

            for (auto & reference_variant_idx: graph[vertex_idx].reference_variant_indices) {

                assert(reference_variant_idx != graph[vertex_idx].variant_allele_idx.first);
                nested_variant_cluster_dependency_it.first->second.push_back(reference_variant_idx);
            }

            sort(nested_variant_cluster_dependency_it.first->second.begin(), nested_variant_cluster_dependency_it.first->second.end(), greater<ushort>());
        }
    }

    return variant_cluster_haplotypes;
}

void VariantClusterGraph::updateVariantPathIndices(vector<pair<ushort, vector<bool> > > * variant_path_indices, unordered_map<pair<ushort, ushort>, pair<uint, uint>, boost::hash<pair<ushort, ushort> > > * running_variants, const uint num_nucleotides, const ushort num_best_paths, const ushort best_paths_idx) {

    auto running_variants_it = running_variants->begin();

    while (running_variants_it != running_variants->end()) {

        assert(running_variants_it->second.first < running_variants_it->second.second);

        if (running_variants_it->second.second <= num_nucleotides) {

            running_variants_it = running_variants->erase(running_variants_it);
            continue;
        }

        if (running_variants_it->second.first <= num_nucleotides) {

            if (variant_path_indices->empty()) {

                variant_path_indices->emplace_back(running_variants_it->first.first, vector<bool>(num_best_paths, false));
                variant_path_indices->back().second.at(best_paths_idx) = true;

            } else {

                bool new_variant_idx = true;

                for (auto & variant_path_idx: *variant_path_indices) {

                    if (variant_path_idx.first == running_variants_it->first.first) {

                        assert(variant_path_idx.second.size() == num_best_paths);
                        variant_path_idx.second.at(best_paths_idx) = true;

                        new_variant_idx = false;
                        break;
                    }
                }

                if (new_variant_idx) {

                    variant_path_indices->emplace_back(running_variants_it->first.first, vector<bool>(num_best_paths, false));
                    variant_path_indices->back().second.at(best_paths_idx) = true;
                }
            }
        } 

        running_variants_it++;
    } 
}

bool VariantClusterGraphCompare(VariantClusterGraph * first_variant_cluster_graph, VariantClusterGraph * second_variant_cluster_graph) { 

    assert(first_variant_cluster_graph);
    assert(second_variant_cluster_graph);

    return (first_variant_cluster_graph->getInfo().size() > second_variant_cluster_graph->getInfo().size());
}


