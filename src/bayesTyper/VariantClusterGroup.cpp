
/*
VariantClusterGroup.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


VariantClusterGroup::VariantClusterGroup(const vector<VariantCluster *> & variant_clusters, const vector<VariantClusterGraph *> & variant_cluster_graphs, const unordered_map<uint, uint> & variant_cluster_depedencies) {

	assert(!(variant_clusters.empty()));
	assert(variant_clusters.size() < Utils::uint_overflow);
	assert(variant_clusters.size() == variant_cluster_graphs.size());
    assert(variant_cluster_depedencies.size() < variant_clusters.size());

    chrom_name = variant_clusters.front()->chrom_name;
    chrom_class = variant_clusters.front()->chrom_class;

	start_position = Utils::uint_overflow;
	end_position = 0;

 	num_variants = 0;

	is_snv = false;

	if ((variant_clusters.size() == 1) and (variant_clusters.front()->variants.size() == 1) and (variant_clusters.front()->variants.begin()->second.type == VariantCluster::VariantType::SNP)) {

		is_snv = true;
	}

	unordered_map<uint, uint> variant_cluster_idx_to_vertex_id;
	vertices.reserve(variant_clusters.size());

	for (uint i = 0; i < variant_clusters.size(); i++) {

		if (variant_cluster_depedencies.count(variant_clusters.at(i)->cluster_idx) < 1) {

			source_vertices.push_back(vertices.size());
		}

		assert(variant_cluster_idx_to_vertex_id.emplace(variant_clusters.at(i)->cluster_idx, vertices.size()).second);

		vertices.emplace_back(variant_clusters.at(i)->cluster_idx);
		vertices.back().graph = variant_cluster_graphs.at(i);

		assert(vertices.back().graph != nullptr);
		assert(vertices.back().genotyper == nullptr);

		assert(variant_clusters.at(i)->contained_clusters.empty());

		assert(chrom_name == variant_clusters.at(i)->chrom_name);
		assert(chrom_class == variant_clusters.at(i)->chrom_class);

		assert(variant_clusters.at(i)->left_flank == variant_clusters.at(i)->variants.begin()->first);
		assert(variant_clusters.at(i)->right_flank >= variant_clusters.at(i)->variants.rbegin()->first);

		start_position = min(start_position, variant_clusters.at(i)->left_flank + 1);
		end_position = max(end_position, variant_clusters.at(i)->right_flank + 1);

		num_variants += variant_cluster_graphs.at(i)->getInfo().size();
	}

	assert(start_position <= end_position);

	out_edges = vector<vector<uint> >(vertices.size());

	for (auto & variant_cluster_depedency: variant_cluster_depedencies) {

		auto vertex_idx_source = variant_cluster_idx_to_vertex_id.at(variant_cluster_depedency.second);
		auto vertex_idx_target = variant_cluster_idx_to_vertex_id.at(variant_cluster_depedency.first);

		assert(find(out_edges.at(vertex_idx_source).begin(), out_edges.at(vertex_idx_source).end(), vertex_idx_target) == out_edges.at(vertex_idx_source).end()); 

		out_edges.at(vertex_idx_source).push_back(vertex_idx_target);
	}

	assert(source_vertices.size() <= vertices.size());
}

VariantClusterGroup::~VariantClusterGroup() {

	for (auto & vertex: vertices) {

		delete vertex.graph;
		delete vertex.genotyper;
	}
}

uint VariantClusterGroup::numberOfVariants() const {

	return num_variants;
}

uint VariantClusterGroup::numberOfVariantClusters() const {

	return vertices.size();
}

uint VariantClusterGroup::numberOfVariantClusterGroupTrees() const {

	return source_vertices.size();
}

void VariantClusterGroup::findSamplePaths(KmerBloom<Utils::kmer_size> * sample_kmer_bloom, const uint prng_seed, const ushort max_sample_haplotype_candidates) {

	for (auto & vertex: vertices) {

		vertex.graph->findSamplePaths(sample_kmer_bloom, prng_seed + vertex.variant_cluster_idx, max_sample_haplotype_candidates);
	}
}

void VariantClusterGroup::countPathKmers(unordered_set<bitset<Utils::kmer_size * 2> > * path_kmers) {

	for (auto & vertex: vertices) {

		vertex.graph->countPathKmers(path_kmers);
	}
}

void VariantClusterGroup::classifyPathKmers(KmerCountsHash * kmer_hash, KmerBloom<Utils::kmer_size> * multigroup_kmer_bloom) {

	for (auto & vertex: vertices) {

		vertex.graph->classifyPathKmers(kmer_hash, multigroup_kmer_bloom);
	}
}

bool VariantClusterGroup::isAutosomalSimpleSNV() const {

	if (!is_snv or (chrom_class != Utils::ChromClass::Autosomal)) {

		return false;
	}

	assert(numberOfVariants() == 1);
	assert(numberOfVariantClusters() == 1);
	assert(numberOfVariantClusterGroupTrees() == 1);

	return vertices.front().graph->isSimpleCluster();
}

void VariantClusterGroup::initGenotyper(KmerCountsHash * kmer_hash, const vector<Sample> & samples, const uint prng_seed, const uchar num_genomic_rate_gc_bias_bins, const float kmer_subsampling_rate, const uint max_haplotype_variant_kmers) {

	for (auto & vertex: vertices) {

		assert(!(vertex.genotyper));

		vertex.genotyper = new VariantClusterGenotyper(vertex.graph, kmer_hash, samples, prng_seed + vertex.variant_cluster_idx, num_genomic_rate_gc_bias_bins);
		vertex.genotyper->reset(kmer_subsampling_rate, max_haplotype_variant_kmers);
	}
}

void VariantClusterGroup::resetGenotyper(const float kmer_subsampling_rate, const uint max_haplotype_variant_kmers) {

	for (auto & vertex: vertices) {

		assert(vertex.genotyper);
		vertex.genotyper->reset(kmer_subsampling_rate, max_haplotype_variant_kmers);
	}
}

void VariantClusterGroup::resetGroup() {

	for (auto & vertex: vertices) {

		assert(vertex.genotyper);
		
		delete vertex.genotyper;
		vertex.genotyper = nullptr;
	}
}

void VariantClusterGroup::shuffleBranchOrdering(const uint prng_seed) {

    mt19937 prng = mt19937(prng_seed);

	shuffle(source_vertices.begin(), source_vertices.end(), prng);

	for (auto & vertex_out_edges: out_edges) {
		
		shuffle(vertex_out_edges.begin(), vertex_out_edges.end(), prng);		
	}
}

void VariantClusterGroup::estimateGenotypes(const CountDistribution & count_distribution, const ChromosomePloidy & chrom_ploidy, const bool collect_samples) {

	vector<VariantClusterHaplotypes::NestedVariantClusterInfo> nested_variant_cluster_info;
	nested_variant_cluster_info.reserve(chrom_ploidy.getPloidy(chrom_class).size());
	
	for (auto & ploidy: chrom_ploidy.getPloidy(chrom_class)) {

		nested_variant_cluster_info.emplace_back(ploidy);
	}

	for (auto & source_vertex: source_vertices) {

		runGibbsSample(source_vertex, count_distribution, nested_variant_cluster_info, collect_samples);
	}
}

void VariantClusterGroup::runGibbsSample(const uint vertex_idx, const CountDistribution & count_distribution, const vector<VariantClusterHaplotypes::NestedVariantClusterInfo> & nested_variant_cluster_info, const bool collect_samples) {

	assert(vertices.at(vertex_idx).genotyper);

	vertices.at(vertex_idx).genotyper->sampleDiplotypes(count_distribution, nested_variant_cluster_info, collect_samples);
	vertices.at(vertex_idx).genotyper->sampleHaplotypeFrequencies();
		
	for (auto & target_vertex_idx: out_edges.at(vertex_idx)) {

		auto target_nested_variant_cluster_info = nested_variant_cluster_info;
		vertices.at(vertex_idx).genotyper->updateNestedVariantClusterInfo(&target_nested_variant_cluster_info, vertices.at(target_vertex_idx).variant_cluster_idx);

		runGibbsSample(target_vertex_idx, count_distribution, target_nested_variant_cluster_info, collect_samples);
	}
}

void VariantClusterGroup::getNoiseCounts(CountAllocation * noise_counts, const CountDistribution & count_distribution) {

	assert(vertices.size() == 1);
	assert(source_vertices.size() == 1);

	assert(vertices.front().genotyper);
	vertices.front().genotyper->getNoiseCounts(noise_counts, count_distribution);	
}

void VariantClusterGroup::collectGenotypes(vector<Genotypes*> * variant_genotypes, const ChromosomePloidy & chrom_ploidy, const Filters & filters) {

	for (auto & vertex: vertices) {	

		assert(vertex.genotyper);
		auto vertex_genotypes = vertex.genotyper->getGenotypes(chrom_name, chrom_ploidy.getPloidy(chrom_class), filters);

		for (auto & genotypes: vertex_genotypes) {

			genotypes->variant_cluster_group_size = vertices.size();
			genotypes->variant_cluster_group_region = chrom_name + ":" + to_string(start_position) + "-" + to_string(end_position);

			variant_genotypes->push_back(genotypes);
		}
	}
}

bool VariantClusterGroupCompare(VariantClusterGroup * first_variant_cluster_group, VariantClusterGroup * second_variant_cluster_group) { 

    assert(first_variant_cluster_group);
    assert(second_variant_cluster_group);

    return (first_variant_cluster_group->numberOfVariants() > second_variant_cluster_group->numberOfVariants());
}


