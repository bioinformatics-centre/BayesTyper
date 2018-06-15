
/*
CountDistribution.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#include <fstream>
#include <math.h>
#include <algorithm>

#include "boost/math/special_functions/factorials.hpp"
#include "boost/math/special_functions/gamma.hpp"

#include "CountDistribution.hpp"
#include "Sample.hpp"
#include "OptionsContainer.hpp"
#include "Combinator.hpp"


CountDistribution::CountDistribution(const vector<Sample> & samples_in, const OptionsContainer & options_container) : samples(samples_in), num_genomic_rate_gc_bias_bins(1), noise_rate_priors(samples.size(), options_container.getValue<pair<double,double> >("noise-rate-prior")) {

	prng = mt19937(options_container.getValue<uint>("random-seed"));
	
	genomic_count_distributions = vector<vector<NegativeBinomialDistribution> >(samples.size(), vector<NegativeBinomialDistribution>(num_genomic_rate_gc_bias_bins));

	noise_rates = vector<double>(samples.size(), 0);
	genomic_count_log_pmf_cache = vector<vector<vector<vector<double> > > >(samples.size(), vector<vector<vector<double> > >(num_genomic_rate_gc_bias_bins, vector<vector<double> >(max_multiplicity + 1, vector<double>(max_kmer_count + 1))));
	noise_count_log_pmf_cache = vector<vector<double> >(samples.size(), vector<double>(max_kmer_count + 1));

	resetNoiseRates();
	updateGenomicCache();
}


void CountDistribution::setGenomicCountDistributions(const vector<vector<KmerStats> > & intercluster_diploid_kmer_stats, const string & output_prefix) {

	assert(intercluster_diploid_kmer_stats.size() == samples.size());
	assert(intercluster_diploid_kmer_stats.size() == genomic_count_distributions.size());

	assert(!(intercluster_diploid_kmer_stats.empty()));
	assert(intercluster_diploid_kmer_stats.front().size() == 1);

    cout << "[" << Utils::getLocalTime() << "] Estimating genomic haploid kmer count distribtion(s) from " << intercluster_diploid_kmer_stats.front().front().getCount() << " diploid parameter kmers ...\n" << endl;

    ofstream genomic_outfile(output_prefix + ".txt");
    assert(genomic_outfile.is_open());

    genomic_outfile << "Sample\tMean\tVariance" << endl;

    for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

    	assert(num_genomic_rate_gc_bias_bins == 1);
		assert(genomic_count_distributions.at(sample_idx).size() == num_genomic_rate_gc_bias_bins);

        for (ushort bias_idx = 0; bias_idx < num_genomic_rate_gc_bias_bins; bias_idx++) {

        	assert(intercluster_diploid_kmer_stats.front().front().getCount() == intercluster_diploid_kmer_stats.at(sample_idx).at(bias_idx).getCount());

        	auto sample_mean = intercluster_diploid_kmer_stats.at(sample_idx).at(bias_idx).getMean();
        	assert(sample_mean.second);

        	auto sample_var = intercluster_diploid_kmer_stats.at(sample_idx).at(bias_idx).getVariance();
        	assert(sample_var.second);

        	auto ng_parameters = NegativeBinomialDistribution::momentsToParameters(sample_mean.first, sample_var.first);
        	ng_parameters.second /= 2;

            genomic_count_distributions.at(sample_idx).at(bias_idx).setParameters(ng_parameters);

            auto ng_mean = genomic_count_distributions.at(sample_idx).at(bias_idx).mean();
            auto ng_var = genomic_count_distributions.at(sample_idx).at(bias_idx).var();

            cout << "[" << Utils::getLocalTime() << "] Estimated fixed negative binomial distribution for sample " << samples.at(sample_idx).name << " with mean " << ng_mean << " and variance " << ng_var << endl;

	        genomic_outfile << samples.at(sample_idx).name << "\t" << ng_mean << "\t" << ng_var << endl;
        }
    }

    genomic_outfile.close();
    updateGenomicCache();

	cout << "\n[" << Utils::getLocalTime() << "] Wrote parameters to " << output_prefix << ".txt" << endl;
}

const vector<vector<NegativeBinomialDistribution> > & CountDistribution::getGenomicCountDistributions() const {

	return genomic_count_distributions;
}

const vector<double> & CountDistribution::getNoiseRates() const {

	return noise_rates;
}

void CountDistribution::setNoiseRates(const vector<double> & noise_rates_in) {

	assert(noise_rates_in.size() == samples.size());
	assert(noise_rates_in.size() == noise_rates.size());
	
	noise_rates = noise_rates_in;
	
	updateNoiseCache();
}

void CountDistribution::resetNoiseRates() {

	for (uint sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

		noise_rates.at(sample_idx) = sampleGamma(noise_rate_priors.at(sample_idx).first, noise_rate_priors.at(sample_idx).second);
	}

	updateNoiseCache();
}

void CountDistribution::sampleNoiseParameters(const CountAllocation & noise_counts) {

	assert(noise_counts.getCounts().size() == samples.size());

	for (uint sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

		assert(noise_counts.getCounts().at(sample_idx).size() == (Utils::uchar_overflow + 1));
		auto noise_count_suff_stats = calcCountSuffStats(noise_counts.getCounts().at(sample_idx));

		noise_rates.at(sample_idx) = sampleGamma(noise_rate_priors.at(sample_idx).first + noise_count_suff_stats.second, noise_rate_priors.at(sample_idx).second/(noise_count_suff_stats.first * noise_rate_priors.at(sample_idx).second + 1));
	}

	updateNoiseCache();
}

pair<ulong, ulong> CountDistribution::calcCountSuffStats(const vector<ulong> & counts) const {

	ulong num_observations = 0;
	ulong count_sum = 0;

	for (uint i = 0; i < counts.size(); i++) {

		num_observations += counts.at(i);
		count_sum += (i * counts.at(i));
	}

	return make_pair(num_observations, count_sum);
}

double CountDistribution::sampleGamma(const double shape, const double scale) {

	assert(shape > 0);
	assert(scale > 0);

	gamma_dist.param(gamma_distribution<>::param_type(shape, scale));
	double sample = gamma_dist(prng);

	assert(sample > 0);

	return sample;
}

void CountDistribution::updateGenomicCache() {

	assert(genomic_count_log_pmf_cache.size() == samples.size());

	for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

		assert(genomic_count_log_pmf_cache.at(sample_idx).size() == num_genomic_rate_gc_bias_bins);

    	for (ushort bias_idx = 0; bias_idx < num_genomic_rate_gc_bias_bins; bias_idx++) {

			assert(genomic_count_log_pmf_cache.at(sample_idx).at(bias_idx).size() == static_cast<uint>(max_multiplicity + 1));

			for (ushort kmer_multiplicity = 0; kmer_multiplicity <= max_multiplicity; kmer_multiplicity++) {

				assert(genomic_count_log_pmf_cache.at(sample_idx).at(bias_idx).at(kmer_multiplicity).size() == static_cast<uint>(max_kmer_count + 1));

				for (ushort kmer_count = 0; kmer_count <= max_kmer_count; kmer_count++) {

					genomic_count_log_pmf_cache.at(sample_idx).at(bias_idx).at(kmer_multiplicity).at(kmer_count) = genomicCountLogPmf(sample_idx, bias_idx, kmer_multiplicity, kmer_count);
				}
			}
		}
	}
}

void CountDistribution::updateNoiseCache() {

	assert(noise_count_log_pmf_cache.size() == samples.size());

	for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

		assert(noise_count_log_pmf_cache.at(sample_idx).size() == static_cast<uint>(max_kmer_count + 1));

		for (ushort kmer_count = 0; kmer_count <= max_kmer_count; kmer_count++) {

			noise_count_log_pmf_cache.at(sample_idx).at(kmer_count) = noiseCountLogPmf(sample_idx, kmer_count);
		}
	}
}

double CountDistribution::calcCountLogProb(const ushort sample_idx, const uchar bias_idx, const uchar kmer_multiplicity, const uchar kmer_count) const {

	if (kmer_multiplicity == 0) {

		return noise_count_log_pmf_cache.at(sample_idx).at(kmer_count);		

	} else {

		return genomic_count_log_pmf_cache.at(sample_idx).at(bias_idx).at(kmer_multiplicity).at(kmer_count);
	}
}

double CountDistribution::genomicCountLogPmf(const ushort sample_idx, const uchar bias_idx, const uchar kmer_multiplicity, const uchar kmer_count) const {

	assert(kmer_count <= max_kmer_count);

	if (kmer_multiplicity == 0) {

		if (kmer_count == 0) {

			return 0;
		
		} else {

			return -std::numeric_limits<double>::infinity();
		} 
	} 

	double genomic_count_log_pmf = genomic_count_distributions.at(sample_idx).at(bias_idx).logPmf(kmer_count, kmer_multiplicity);

	if (kmer_count == max_kmer_count) {

		uint kmer_count_limit = kmer_count;
		double prev_genomic_count_log_pmf = 0;

		do {

			kmer_count_limit++;
			prev_genomic_count_log_pmf = genomic_count_log_pmf;

			genomic_count_log_pmf = Utils::logAddition(genomic_count_log_pmf, genomic_count_distributions.at(sample_idx).at(bias_idx).logPmf(kmer_count_limit, kmer_multiplicity));

			if (genomic_count_log_pmf > 0) {

				assert(isfinite(genomic_count_log_pmf));
				genomic_count_log_pmf = 0;

				break;
			}

		} while (!(Utils::doubleCompare(prev_genomic_count_log_pmf, genomic_count_log_pmf)));
	}

	assert(isfinite(genomic_count_log_pmf));
	assert(genomic_count_log_pmf <= 0);

	return genomic_count_log_pmf;
}

double CountDistribution::noiseCountLogPmf(const ushort sample_idx, const uchar kmer_count) const {

	assert(kmer_count <= max_kmer_count);

	double noise_count_log_pmf = poissonLogProb(kmer_count, noise_rates.at(sample_idx));

	if (kmer_count == max_kmer_count) {

		uint kmer_count_limit = kmer_count;
		double prev_noise_count_log_pmf = 0;

		do {

			kmer_count_limit++;
			prev_noise_count_log_pmf = noise_count_log_pmf;

			noise_count_log_pmf = Utils::logAddition(noise_count_log_pmf, poissonLogProb(kmer_count_limit, noise_rates.at(sample_idx)));

			if (noise_count_log_pmf > 0) {

				assert(isfinite(noise_count_log_pmf));
				noise_count_log_pmf = 0;

				break;
			}

		} while (!(Utils::doubleCompare(prev_noise_count_log_pmf, noise_count_log_pmf)));
	} 

	assert(isfinite(noise_count_log_pmf));
	assert(noise_count_log_pmf <= 0);

	return noise_count_log_pmf;
}

double CountDistribution::poissonLogProb(const uint value, const double rate) const {

	return value * log(rate) - rate - boost::math::lgamma(value + 1);
}


