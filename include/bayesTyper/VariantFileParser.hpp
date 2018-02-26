
/*
VariantFileParser.hpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#ifndef __bayesTyper__VariantFileParser_hpp
#define __bayesTyper__VariantFileParser_hpp

#include <string>
#include <unordered_map>
#include <list>
#include <mutex>
#include <condition_variable>
#include <memory>
#include <thread>
#include <set>

#include "boost/graph/adjacency_list.hpp"
#include "ProducerConsumerQueue.hpp"

#include "Utils.hpp"
#include "VariantCluster.hpp"
#include "VariantClusterGraph.hpp"
#include "VariantClusterGroup.hpp"

using namespace std;


class VariantFileParser {
  
	public:

		VariantFileParser(const string, const uint, const double, const ushort);
		~VariantFileParser();

		void addDecoys(const string, const uchar);

		template <uchar kmer_size>
		void readVariantFile(const string, vector<VariantClusterGraph*> *, vector<VariantClusterGroup*> *, const uint);

		struct InterClusterRegion {

			string::const_iterator start;
			string::const_iterator end;
			Utils::ChromosomeClass chromosome_class;

			InterClusterRegion() {}
			InterClusterRegion(string::const_iterator start_in, string::const_iterator end_in, Utils::ChromosomeClass chromosome_class_in) : start(start_in), end(end_in), chromosome_class(chromosome_class_in) {}
		};
		
		const vector<InterClusterRegion> & getInterclusterRegions();
		
		ulong getNumberOfVariants();
		ushort getMaxAlternativeAlleles();
		ulong getVariableRegionLength();

		static string simplifyChromosomeId(const string & chromosome);
		static Utils::ChromosomeClass classifyGenomeChromosome(const string & chromosome);

	private:

		enum class AlleleCount : uchar {Total = 0, Excluded_genome, Excluded_match, Excluded_end, Excluded_length, ALLELE_COUNT_SIZE};

		const uint max_allele_length;
		const double copy_number_variant_threshold;
		const ushort num_threads;

		ulong number_of_variants;
		ushort max_alternative_alleles;
		ulong variable_region_length;

		unordered_map<string, string *> genome_sequences;
		unordered_map<string, string *> decoy_sequences;

		vector<InterClusterRegion> intercluster_regions;

	    map<uint, ulong> variant_cluster_size_stats;
	    map<uint, ulong> variant_cluster_group_size_stats;

		ulong parseFasta(const string, unordered_map<string, string *> *);

		void addSequenceToInterclusterRegions(const string, string::const_iterator, string::const_iterator, const uchar);

		template <typename FileType, uchar kmer_size>
		pair<vector<uint>, vector<uint> > parseVariants(FileType *, ProducerConsumerQueue<vector<unordered_map<uint, VariantCluster*> * > * > *);

		void rightTrimAllele(string *, string *);
		void parseAllele(const string &, const string &, const ushort, VariantCluster::Variant *);
		Utils::VariantType classifyAllele(const int, const int);		

		template <uchar kmer_size>
		uint copyNumberVariantLength(const string &, const uint, const string &, const uint);

		template <uchar kmer_size>
		void clusterVariants(VariantCluster::Variant &, const uint, const set<uint>, const string &, const string &, uint *, map<uint, VariantCluster*> *, unordered_map<uint, VariantCluster*> *, list<unordered_set<uint> > *);
		
		void processVariantClusters(ProducerConsumerQueue<vector<unordered_map<uint, VariantCluster*> * > * > *, vector<unordered_map<uint, VariantCluster*> * > **, uint *, unordered_map<uint, VariantCluster*> **, list<unordered_set<uint> > *, map<uint, VariantCluster*> *);
		void mergeVariantClusters(unordered_map<uint, VariantCluster*> *, list<unordered_set<uint> > &);
		
		template <uchar kmer_size>
		void processVariantClusterGroupsCallback(vector<VariantClusterGraph*> *, vector<VariantClusterGroup*> *, ProducerConsumerQueue<vector<unordered_map<uint, VariantCluster*> * > * > *, mutex *);

		unordered_map<uint, uint> getVariantClusterGroupDependencies(unordered_map<uint, VariantCluster*> *);
	
		Utils::ChromosomeClass classifyChromosome(const string & chromosome);
};

bool InterclusterRegionsCompare(const VariantFileParser::InterClusterRegion &, const VariantFileParser::InterClusterRegion &);


#endif
