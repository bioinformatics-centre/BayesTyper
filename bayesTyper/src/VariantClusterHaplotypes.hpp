
/*
VariantClusterHaplotypes.hpp - This file is part of BayesTyper (v0.9)


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


#ifndef __bayesTyper__VariantClusterHaplotypes_hpp
#define __bayesTyper__VariantClusterHaplotypes_hpp

#include <vector>
#include <unordered_map>
#include <random>

#include "../Eigen/Dense"

#include "boost/functional/hash.hpp"

#include "Utils.hpp"
#include "KmerCounts.hpp"
#include "KmerCoverage.hpp"
#include "Sample.hpp"

typedef unordered_map<pair<ushort,ushort>, vector<double>, boost::hash<pair<ushort, ushort> > > SampleFrequencyContainer;

class VariantClusterHaplotypes {

	public:

		VariantClusterHaplotypes() {};
		~VariantClusterHaplotypes();

		vector<unordered_map<ushort, ushort> > variants;
		vector<unordered_set<uint> > nested_variant_cluster_indices;

		Eigen::MatrixXuchar kmer_haplotype_multiplicities;
		vector<KmerCounts *> kmers;

		vector<uint> unique_kmer_indices;
		vector<uint> multicluster_kmer_indices;
		vector<uint> multicluster_kmer_indices_subset;

		vector<unordered_map<ushort, vector<uint> > > haplotype_variant_kmer_indices;		
		vector<vector<uint> > haplotype_multicluster_kmer_indices;
		
		vector<vector<bool> > redundant_multicluster_haplotypes;

		bool empty();

		uchar getDiplotypeInterclusterKmerMultiplicity(const uint, const pair<ushort,ushort> &, const Utils::Sex);
		uchar getDiplotypeMultilusterKmerMultiplicity(const uint, const pair<ushort,ushort> &, const pair<ushort,ushort> &, const ushort, const Utils::Sex);

 		void resetMulticlusterKmers(const uint, mt19937 *, const ushort);

		bool hasUpdatedMulticlusterMultiplicity(const ushort);
		void updateMulticlusterMultiplicityIndices(const uint, const ushort);		

		bool isMulticlusterMultiplicityConstant(const pair<ushort,ushort> &, const pair<ushort,ushort> &);

		void removeDiplotypeMulticlusterMultiplicities(const pair<ushort,ushort> &, const ushort, const uint);
		void addDiplotypeMulticlusterMultiplicities(const pair<ushort,ushort> &, const pair<ushort,ushort> &, const ushort, const uint);
		
		KmerCoverage getAlleleKmerCoverage(const pair<ushort, ushort> &, const ushort, const vector<pair<ushort, ushort> > &, const ushort);

	private:

		uchar getDiplotypeKmerMultiplicity(const uint, const pair<ushort,ushort> &);

		void removeDiplotypeMulticlusterMultiplicity(const uint, const pair<ushort,ushort> &, const ushort, const uint);
		void addDiplotypeMulticlusterMultiplicity(const uint, const pair<ushort,ushort> &, const pair<ushort,ushort> &, const ushort, const uint);

		KmerCoverage getAllelePairKmerCoverage(const ushort, const ushort, const ushort, const ushort);
		void updateAlleleKmerCoverage(KmerCoverage *, KmerCounts *, const ushort);
};

#endif