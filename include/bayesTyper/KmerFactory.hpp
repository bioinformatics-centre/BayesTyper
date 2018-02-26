
/*
KmerFactory.hpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#ifndef __bayesTyper__KmerFactory_hpp
#define __bayesTyper__KmerFactory_hpp

#include <map>
#include <unordered_set>
#include <vector>
#include <string>
#include <random>

#include "boost/graph/adjacency_list.hpp"

#include "Utils.hpp"
#include "KmerHash.hpp"
#include "VariantFileParser.hpp"
#include "VariantClusterGroup.hpp"
#include "Sample.hpp"
#include "OptionsContainer.hpp"
#include "NegativeBinomialDistribution.hpp"

using namespace std;

class KmerFactory {

	private:

		const uint prng_seed;
		const string vcf_file;
		const string genome_file;
		const string output_prefix;
		const string decoy_file;
		const ushort num_threads;
		const uint max_allele_length;
		const double copy_number_variant_threshold;
		const ushort max_sample_haplotype_candidates;
		const uchar num_genomic_rate_gc_bias_bins;

		ulong number_of_variants;
		ushort max_alternative_alleles;

		vector<vector<NegativeBinomialDistribution> > genomic_count_distributions;

	public:
		
		KmerFactory(const OptionsContainer &);

		template <uchar kmer_size>
		KmerHash * initKmerHash(vector<VariantClusterGroup*> *, const vector<Sample> &);

		void estimateGenomicCountDistributions(const vector<vector<vector<ulong> > > &, const vector<Sample> &);

		ulong numberOfVariants();
		ushort maxAlternativeAlleles();

		const vector<vector<NegativeBinomialDistribution> > & genomicCountDistributions();
};


#endif