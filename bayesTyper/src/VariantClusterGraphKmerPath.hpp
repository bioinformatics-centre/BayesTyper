
/*
VariantClusterGraphKmerPath.hpp - This file is part of BayesTyper (v0.9)


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


#ifndef __bayesTyper__VariantClusterGraphKmerPath_hpp
#define __bayesTyper__VariantClusterGraphKmerPath_hpp

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <memory>
#include <algorithm>

#include "boost/functional/hash.hpp"

#include "Utils.hpp"
#include "Kmer.hpp"
#include "KmerCounts.hpp"
#include "LinearMap.hpp"
#include "VariantKmerStats.hpp"
#include "VariantClusterGraphVertex.hpp"
#include "KmerInfo.hpp"


using namespace std;

template<uchar kmer_size>
class VariantClusterGraphKmerPath {
	
	private:
		
		typedef vector<pair<bitset<kmer_size * 2>, KmerInfo> > KmerMultiplicities;
	
		const bool is_dummy;

		double calculateVariantScore(pair<uint, uint> &);
		pair<double, uint> calculatePathScore(const unordered_set<const VariantClusterGraphVertex *> &);

	public:

		VariantClusterGraphKmerPath(const uint, const bool);
		virtual ~VariantClusterGraphKmerPath() {};

		bool operator == (const VariantClusterGraphKmerPath<kmer_size> & rhs) const;

		void addVertex(VariantClusterGraphVertex &, const uint);
		virtual void newKmer(KmerCounts *, const ushort, const bool);	

		const vector<KmerMultiplicities *> & getKmerVertexMultiplicities();
		const vector<const VariantClusterGraphVertex *> & getPath();

		const unordered_map<pair<ushort, ushort>, const VariantClusterGraphVertex *, boost::hash<pair<ushort, ushort> > > & getVariants();
		uint getNumberOfVariantKmers(const ushort);

		uint getNumberOfReferenceAlleles();

		KmerMultiplicities * newKmerVertex(const uint);
		void addKmerVertex(KmerMultiplicities *);

		void incrementNucleotide();		

		bool isDummy();

		pair<double, uint> getScore();
		pair<double, uint> getScore(const unordered_set<const VariantClusterGraphVertex *> &);

		void addCoveredVertices(unordered_set<const VariantClusterGraphVertex *> *);

		KmerPair<kmer_size> kmer_pair;

	protected:

		uint num_nucleotides;
		uint num_reference_alleles;

		vector<KmerMultiplicities *> kmers;
		vector<const VariantClusterGraphVertex *> path;

		unordered_map<pair<ushort, ushort>, const VariantClusterGraphVertex *, boost::hash<pair<ushort, ushort> > > variants;

		struct VariantScore {

			uint last_nucleotide;
			pair<uint, uint> score_counts;

			VariantScore(const uint last_nucleotide_in) : last_nucleotide(last_nucleotide_in) , score_counts(make_pair(0, 0)) {}
		};

		unordered_map<ushort, VariantScore> running_variants;
		unordered_map<ushort, VariantScore> processed_variants;

		void updateScore(KmerCounts *, const ushort);

		virtual void addNestedVariantClusterIndices(VariantClusterGraphVertex &);

};
 

template<uchar kmer_size>
class VariantClusterGraphFullKmerPath : public VariantClusterGraphKmerPath<kmer_size> {

	public:

		VariantClusterGraphFullKmerPath(const uint, const bool);

		void newKmer(KmerCounts *, const ushort, const bool);	
		const unordered_set<uint> & getNestedVariantClusterIndices();		
		uint getNumberOfMultiClusterKmers();		

	private:

		uint num_multicluster_kmers;
		unordered_set<uint> nested_variant_cluster_indices;

		void addNestedVariantClusterIndices(VariantClusterGraphVertex &);
};


#include "VariantClusterGraphKmerPath.tpp"


#endif