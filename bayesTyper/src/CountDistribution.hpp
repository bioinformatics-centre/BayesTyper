
/*
CountDistribution.hpp - This file is part of BayesTyper (v0.9)


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

#include "../Eigen/Dense"

#include "Utils.hpp"
#include "Combinator.hpp"
#include "CountAllocation.hpp"
#include "DiscreteSampler.hpp"
#include "Sample.hpp"
#include "OptionsContainer.hpp"
#include "NegativeBinomialDistribution.hpp"


class CountDistribution {

	public:

		CountDistribution(const ushort, const ushort, const vector<NegativeBinomialDistribution> &, const OptionsContainer &);

		const vector<NegativeBinomialDistribution> & genomicCountDistributions() const;
		const vector<vector<double> > & noiseZeroInflations() const;
		const vector<vector<double> > & noiseRates() const;

		const vector<vector<vector<double> > > & noiseZeroInflationsSamples() const;
		const vector<vector<vector<double> > > & noiseRatesSamples() const;

		void setGenomicCountDistributions(const vector<NegativeBinomialDistribution> &);
		void setNoiseZeroInflations(const vector<vector<double> > &);
		void setNoiseRates(const vector<vector<double> > &);

		void updateMaxMultiplicityCache(const uchar);

		void writeParameterSamples(const vector<Sample> &) const;
		void writeNoiseSamples(const vector<vector<vector<double> > > &, const string &, const ushort, const vector<vector<pair<double,double> > > &, const vector<vector<double> > &, const vector<Sample> &) const;

		double calcCombinedGenotypeNoiseCountLogProb(const ushort, const uchar, const uchar) const;
		double calcCombinedNoiseCountLogProb(const ushort, const uchar) const;

		unique_ptr<LogDiscreteSampler> calcGenotypeCountLogPosterior(const ushort, const uchar, const uchar) const;
		unique_ptr<LogDiscreteSampler> calcNoiseSplitCountLogPosterior(const ushort, const uchar, const uchar) const;

		void sampleParameters(const CountAllocation &);

	private:

		double calcCombinedGenotypeNoiseCountLogProbInternal(const ushort, const uchar, const uchar) const;
		double calcCombinedNoiseCountLogProbInternal(const ushort, const vector<uchar> &, const uchar) const;

		double genomicCountLogPmf(const ushort, const uchar, const uchar) const;
		double genomicCountLogPmfInternal(const ushort, const uchar, const uchar) const;
		double poissonLogPmf(const double, const uchar) const;
		double noiseCountLogPmf(const ushort, const uchar, const uchar) const;

		pair<ulong,ulong> calcCountSuffStats(const unordered_map<uchar,ulong> &) const;

		void sampleNoiseParameters(const vector<vector<unordered_map<uchar,ulong> > > &);

		double sampleGamma(const double, const double);
		double sampleBeta(const double, const double);

		void updateCache();

		const ushort num_samples;
		const ushort num_noise_sources;

		const vector<vector<pair<double,double> > > noise_zero_inflation_priors;
		const vector<vector<pair<double,double> > > noise_rate_priors;

		vector<NegativeBinomialDistribution> genomic_count_distributions;
		vector<vector<double> > noise_zero_inflations;
		vector<vector<double> > noise_rates;

		uchar max_multiplicity_cache;
		uchar max_kmer_count;

		vector<vector<vector<double> > > genomic_count_log_pmf_cache;
		vector<vector<double> > combined_noise_count_log_prob_cache;
		vector<vector<vector<double> > > combined_genotype_noise_count_log_prob_cache;

		mt19937 prng;
		gamma_distribution<> gamma_dist;
		binomial_distribution<> binomial_dist;

		vector<vector<vector<double> > > noise_zero_inflation_samples;
		vector<vector<vector<double> > > noise_rate_samples;
};

#endif
