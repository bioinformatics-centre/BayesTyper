
/*
VariantClusterGroup.cpp - This file is part of BayesTyper (v0.9)


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


#include <unordered_map>
#include <string>
#include <sstream>

#include "boost/graph/adjacency_list.hpp"

#include "VariantClusterGroup.hpp"
#include "Utils.hpp"
#include "VariantClusterGraph.hpp"
#include "VariantClusterGenotyper.hpp"
#include "KmerHash.hpp"
#include "Genotypes.hpp"
#include "ChromosomePloidy.hpp"
#include "Sample.hpp"
#include "VariantInfo.hpp"


VariantClusterGroup::VariantClusterGroup(const vector<VariantCluster *> & variant_clusters, const vector<VariantClusterGraph *> & variant_cluster_graphs, const unordered_map<uint, uint> & variant_cluster_depedencies, const Utils::ChromosomeClass chromosome_class_in, const string & chromosome_id_in) : chromosome_class(chromosome_class_in), chromosome_id(chromosome_id_in) {

	assert(!(variant_clusters.empty()));
	assert(variant_clusters.size() < Utils::uint_overflow);
	assert(variant_clusters.size() == variant_cluster_graphs.size());
    assert(variant_cluster_depedencies.size() < variant_clusters.size());

	start_position = Utils::uint_overflow;
	end_position = 0;

 	number_of_variants = 0;

	is_single_nucleotide_polymorphism = false;

	if ((variant_clusters.size() == 1) and (variant_clusters.front()->variants.size() == 1) and (variant_clusters.front()->variants.begin()->second.type == Utils::VariantType::SNP)) {

		is_single_nucleotide_polymorphism = true;
	}

	map<uint, uint> variant_cluster_idx_to_vertex_id;
	vertices.reserve(variant_clusters.size());

	for (uint i = 0; i < variant_clusters.size(); i++) {

		if (variant_cluster_depedencies.count(variant_clusters.at(i)->cluster_idx) < 1) {

			source_vertices.push_back(vertices.size());
		}

		assert(variant_cluster_idx_to_vertex_id.emplace(variant_clusters.at(i)->cluster_idx, vertices.size()).second);

		vertices.emplace_back(variant_clusters.at(i)->cluster_idx, variant_clusters.at(i)->left_flank + 1, variant_clusters.at(i)->variants.rbegin()->first + 1);

		vertices.back().graph = variant_cluster_graphs.at(i);
		vertices.back().genotyper = nullptr;

		assert(variant_clusters.at(i)->contained_clusters.empty());
		assert(variant_clusters.at(i)->left_flank == variant_clusters.at(i)->variants.begin()->first);

		start_position = min(start_position, vertices.back().start_position);
		end_position = max(end_position, vertices.back().end_position);

		number_of_variants += variant_cluster_graphs.at(i)->getInfo().size();
	}

	for (auto & variant_cluster_depedency: variant_cluster_depedencies) {

		auto out_edges_emplace = out_edges.emplace(variant_cluster_idx_to_vertex_id.at(variant_cluster_depedency.second), vector<uint>());

		assert(find(out_edges_emplace.first->second.begin(), out_edges_emplace.first->second.end(), variant_cluster_idx_to_vertex_id.at(variant_cluster_depedency.first)) == out_edges_emplace.first->second.end()); 
		out_edges_emplace.first->second.push_back(variant_cluster_idx_to_vertex_id.at(variant_cluster_depedency.first));
	}

	complexity = number_of_variants;

	assert(source_vertices.size() <= vertices.size());
}


VariantClusterGroup::~VariantClusterGroup() {

	for (auto & vertex: vertices) {

		delete vertex.graph;
		delete vertex.genotyper;
	}
}


double VariantClusterGroup::getComplexity() const {

	return complexity;
}


uint VariantClusterGroup::numberOfVariants() const {

	return number_of_variants;
}


uint VariantClusterGroup::numberOfVariantClusters() const {

	return vertices.size();
}


uint VariantClusterGroup::numberOfVariantClusterGroupTrees() const {

	return source_vertices.size();
}


bool VariantClusterGroup::isAutosomal() const {

	return (chromosome_class == Utils::ChromosomeClass::Autosomal);
}


bool VariantClusterGroup::isSingleNucleotidePolymorphism() const {

	return is_single_nucleotide_polymorphism;
}


bool VariantClusterGroup::hasAmbiguousNucleotide() const {

	assert(vertices.size() == 1);
	assert(source_vertices.size() == 1);

	assert(!(vertices.front().genotyper));
	return vertices.front().graph->hasAmbiguousNucleotide();
}


bool VariantClusterGroup::hasComplexRegion() const {

	assert(vertices.size() == 1);
	assert(source_vertices.size() == 1);

	assert(!(vertices.front().genotyper));
	return vertices.front().graph->hasComplexRegion();
}


bool VariantClusterGroup::hasRedundantSequence() const {

	assert(vertices.size() == 1);
	assert(source_vertices.size() == 1);

	assert(!(vertices.front().genotyper));
	return vertices.front().graph->hasRedundantSequence();
}


bool VariantClusterGroup::hasNonUniqueKmer() const {

	assert(vertices.size() == 1);
	assert(source_vertices.size() == 1);

	assert(!(vertices.front().genotyper));
	return vertices.front().graph->hasNonUniqueKmer();
}


bool VariantClusterGroup::hasInterclusterKmer() const {

	bool has_intercluster_kmer = false;

	for (auto & vertex: vertices) {

		assert(vertex.genotyper);
		has_intercluster_kmer = vertex.genotyper->hasInterclusterKmer();

		if (has_intercluster_kmer) {

			break;
		}
	}

	return has_intercluster_kmer;
}


bool VariantClusterGroup::hasMulticlusterKmer() const {

	bool has_multicluster_kmer = false;

	for (auto & vertex: vertices) {

		assert(vertex.genotyper);
		has_multicluster_kmer = vertex.genotyper->hasMulticlusterKmer();

		if (has_multicluster_kmer) {

			break;
		}
	}

	return has_multicluster_kmer;
}


bool VariantClusterGroup::hasExcludedKmer() const {

	bool has_excluded_kmer = false;

	for (auto & vertex: vertices) {

		assert(vertex.genotyper);
		has_excluded_kmer = vertex.genotyper->hasExcludedKmer();

		if (has_excluded_kmer) {

			break;
		}
	}

	return has_excluded_kmer;
}


ulong VariantClusterGroup::countSmallmers(Utils::SmallmerSet * smallmer_set) {

	ulong num_unique_smallmers = 0;

	for (auto & vertex: vertices) {

		num_unique_smallmers += vertex.graph->countSmallmers(smallmer_set);
	}

	return num_unique_smallmers;
}


ulong VariantClusterGroup::countKmers(KmerHash * kmer_hash, const uint prng_seed, const ushort num_haplotype_candidates_per_sample, const uint variant_cluster_group_idx) {

	complexity = 0;
	ulong num_unique_kmers = 0;

	for (auto & vertex: vertices) {

		num_unique_kmers += vertex.graph->countKmers(kmer_hash, prng_seed, num_haplotype_candidates_per_sample, variant_cluster_group_idx, &complexity);
	}

	return num_unique_kmers;
}


bool VariantClusterGroup::isInChromosomeRegions(const Regions & chromosome_regions) {

	return chromosome_regions.isIn(chromosome_id, start_position, end_position);
}


void VariantClusterGroup::initialise(KmerHash * kmer_hash, const vector<Sample> & samples, const uint prng_seed, const uchar num_noise_sources, const ushort num_haplotype_candidates_per_sample, const uint max_multicluster_kmers) {

	for (auto & vertex: vertices) {

		if (vertex.genotyper) {

			assert(!(vertex.graph));
			vertex.genotyper->reset(max_multicluster_kmers);
			
		} else {	

			assert(vertex.graph);

			vertex.genotyper = new VariantClusterGenotyper(samples, vertex.variant_cluster_idx);
			vertex.genotyper->initialise(prng_seed, accumulate(chromosome_id.begin(), chromosome_id.end(), start_position + vertex.variant_cluster_idx), kmer_hash, vertex.graph, num_noise_sources, num_haplotype_candidates_per_sample, max_multicluster_kmers);
			
			delete vertex.graph;
			vertex.graph = nullptr;
		}
	}
}


void VariantClusterGroup::shuffleBranches(mt19937 * prng) {

	shuffle(source_vertices.begin(), source_vertices.end(), *prng);

	for (auto & vertex_out_edges: out_edges) {
		
		shuffle(vertex_out_edges.second.begin(), vertex_out_edges.second.end(), *prng);		
	}
}


void VariantClusterGroup::estimateGenotypes(const CountDistribution & count_distribution, const ChromosomePloidy & chromosome_ploidy, const bool collect_samples) {

	for (auto & source_vertex: source_vertices) {

		runGibbsSample(source_vertex, count_distribution, chromosome_ploidy.getPloidy(chromosome_class), collect_samples);
	}
}


void VariantClusterGroup::runGibbsSample(const uint vertex_idx, const CountDistribution & count_distribution, const vector<Utils::Ploidy> & variant_cluster_ploidy, const bool collect_samples) {

	assert(!(vertices.at(vertex_idx).graph));
	assert(vertices.at(vertex_idx).genotyper);

	vertices.at(vertex_idx).genotyper->sampleGenotypes(count_distribution, variant_cluster_ploidy, collect_samples);
	vertices.at(vertex_idx).genotyper->sampleHaplotypeFrequencies();

	auto vertex_out_edges_it = out_edges.find(vertex_idx);

	if (vertex_out_edges_it != out_edges.end()) {
		
		for (auto & next_vertex_idx: vertex_out_edges_it->second) {

			auto next_variant_cluster_ploidy = vertices.at(vertex_idx).genotyper->getVariantClusterPloidy(vertices.at(next_vertex_idx).variant_cluster_idx);
			assert(next_variant_cluster_ploidy.size() == variant_cluster_ploidy.size());

			runGibbsSample(next_vertex_idx, count_distribution, next_variant_cluster_ploidy, collect_samples);
		}
	}
}


void VariantClusterGroup::sampleCountAllocations(const CountDistribution & count_distribution, CountAllocation * count_allocation) {

	assert(vertices.size() == 1);
	assert(source_vertices.size() == 1);

	assert(!(vertices.front().graph));
	assert(vertices.front().genotyper);

	vertices.front().genotyper->sampleCountAllocations(count_distribution);	
	count_allocation->addNoiseCounts(vertices.front().genotyper->getNoiseCounts());
}


void VariantClusterGroup::collectGenotypes(vector<Genotypes*> * variant_genotypes, const ChromosomePloidy & chromosome_ploidy) {

	for (auto & vertex: vertices) {	

		assert(!(vertex.graph));
		assert(vertex.genotyper);

		auto estimated_genotypes = vertex.genotyper->getGenotypes(chromosome_ploidy.getPloidy(chromosome_class));			

		for (auto &git: estimated_genotypes) {

			assert(vertex.start_position <= vertex.end_position);
			assert(start_position <= end_position);

			git->variant_cluster_id = chromosome_id + ":" + to_string(vertex.start_position) + "-" + to_string(vertex.end_position);
			git->variant_cluster_group_id = chromosome_id + ":" + to_string(start_position) + "-" + to_string(end_position);
			git->variant_cluster_group_size = vertices.size();

			variant_genotypes->push_back(git);
			git++;
		}
	}
}


bool VariantClusterGroupCompare(VariantClusterGroup * first_variant_cluster_group, VariantClusterGroup * second_variant_cluster_group) { 

    assert(first_variant_cluster_group);
    assert(second_variant_cluster_group);

    return (first_variant_cluster_group->getComplexity() > second_variant_cluster_group->getComplexity());
}


