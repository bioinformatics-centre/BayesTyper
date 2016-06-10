
/*
VariantFileParser.hpp - This file is part of BayesTyper (v0.9)


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

		VariantFileParser(const string, const uint, const ushort, const ushort, const uint);
		~VariantFileParser();

		template <int kmer_size>
		void addDecoys(const string);

		template <int kmer_size>
		void readVariantFile(const string, vector<VariantClusterGroup*> *);
		
		void writeSizeStatisticsFiles();
		vector<pair<string::const_iterator, string::const_iterator> > & getInterclusterIntervals(const Utils::ChromosomeClass);
		
		ulong getNumberOfVariants();
		ushort getMaxAlternativeAlleles();

		static string simplifyChromosomeId(const string & chromosome);
		static Utils::ChromosomeClass classifyGenomeChromosome(const string & chromosome);	

	private:

		enum class AlleleCount : uchar {Total = 0, Excluded_genome, Excluded_match, Excluded_end, Excluded_length, Excluded_kmers, ALLELE_COUNT_SIZE};

		const uint max_allele_length;
		const uint min_allele_kmers;
		const ushort num_threads;

		mt19937 prng;

		ulong number_of_variants;
		ulong number_of_variant_clusters;
		ushort max_alternative_alleles;

		unordered_map<string, string *> genome_sequences;
		unordered_map<string, string *> decoy_sequences;

		vector<vector<pair<string::const_iterator, string::const_iterator> > > intercluster_intervals;

	    map<uint, ulong> variant_cluster_size_stats;
	    map<uint, ulong> variant_cluster_group_size_stats;

		ulong parseFasta(const string, unordered_map<string, string *> *);

		template <typename FileType, int kmer_size>
		pair<vector<uint>, vector<uint> > parseVariants(FileType *, ProducerConsumerQueue<vector<unordered_map<uint, VariantCluster*> * > * > *);

		uint rightTrimAllele(const string &, string *);

		template <int kmer_size>
		uint numberOfAlleleKmers(const string &, const string &, const uint);

		string reverseComplementSequence(const string &);

		template <int kmer_size>
		uint copyNumberVariantLength(const string &, const uint, const string &, const uint);
		
		template <int kmer_size>
		void parseAllele(string &, string &, const uint, const uint, const ushort, VariantCluster::Variant *, string &);
		
		Utils::VariantType classifyAllele(const int, const int);

		template <int kmer_size>
		void clusterVariants(VariantCluster::Variant &, const uint, const set<uint>, const string &, const string &, uint *, map<uint, VariantCluster*> *, unordered_map<uint, VariantCluster*> *, list<unordered_set<uint> > *);
		
		void processVariantClusters(ProducerConsumerQueue<vector<unordered_map<uint, VariantCluster*> * > * > *, vector<unordered_map<uint, VariantCluster*> * > **, uint *, unordered_map<uint, VariantCluster*> **, list<unordered_set<uint> > *, map<uint, VariantCluster*> *);
		void mergeVariantClusters(unordered_map<uint, VariantCluster*> *, list<unordered_set<uint> > &);
		
		template <int kmer_size>
		void processVariantClusterGroupsCallback(ProducerConsumerQueue<vector<unordered_map<uint, VariantCluster*> * > * > *, mutex *, vector<VariantClusterGroup*> *);

		unordered_map<uint, uint> getVariantClusterGroupDependencies(unordered_map<uint, VariantCluster*> *);
		
		template <int kmer_size>
		VariantClusterGroup * constructVariantClusterGroup(const unordered_map<uint, VariantCluster*> &, const unordered_map<uint, uint> &);
	
		template <int kmer_size>
		void addSequencesToInterclusterIntervals(unordered_map<string, string *> &, unordered_set<string> *);

		void writeSizeStatisticsFile(string, map<uint, ulong> &);
		Utils::ChromosomeClass classifyChromosome(const string & chromosome);
};


#endif
