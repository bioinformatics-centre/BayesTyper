
/*
CountDistribution.hpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#ifndef __bayesTyper__CountDistribution_hpp
#define __bayesTyper__CountDistribution_hpp

#include <list>
#include <unordered_map>
#include <random>
#include <string>
#include <memory>

#include "Utils.hpp"
#include "Combinator.hpp"
#include "CountAllocation.hpp"
#include "DiscreteSampler.hpp"
#include "Sample.hpp"
#include "OptionsContainer.hpp"
#include "NegativeBinomialDistribution.hpp"
#include "KmerStats.hpp"


class CountDistribution {

	public:

		CountDistribution(const vector<Sample> &, const OptionsContainer &);

		const vector<vector<NegativeBinomialDistribution> > & getGenomicCountDistributions() const;
		void setGenomicCountDistributions(const vector<vector<KmerStats> > &, const string &);
		
		const vector<double> & getNoiseRates() const;
		void setNoiseRates(const vector<double> &);
		void resetNoiseRates();

		void sampleNoiseParameters(const CountAllocation &);
		double calcCountLogProb(const ushort, const uchar, const uchar, const uchar) const;

	private:

		double sampleGamma(const double, const double);
		pair<ulong, ulong> calcCountSuffStats(const vector<ulong> &) const;

		void updateGenomicCache();
		void updateNoiseCache();

		double genomicCountLogPmf(const ushort, const uchar, const uchar, const uchar) const;
		double noiseCountLogPmf(const ushort, const uchar) const;

		double poissonLogProb(const uint, const double) const;

		const vector<Sample> samples;	
		const uchar num_genomic_rate_gc_bias_bins;
		
		const uchar max_multiplicity = Utils::uchar_overflow;
		const uchar max_kmer_count = Utils::uchar_overflow;

		vector<vector<NegativeBinomialDistribution> > genomic_count_distributions;

		const vector<pair<double,double> > noise_rate_priors;
		vector<double> noise_rates;

		vector<vector<vector<vector<double> > > > genomic_count_log_pmf_cache;
		vector<vector<double> > noise_count_log_pmf_cache;

		mt19937 prng;
		gamma_distribution<> gamma_dist;
};

#endif
