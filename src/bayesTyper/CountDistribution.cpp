
/*
CountDistribution.cpp - This file is part of BayesTyper (v0.9)


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


CountDistribution::CountDistribution(const ushort num_samples_in, const vector<vector<NegativeBinomialDistribution> > & genomic_count_distributions_in, const OptionsContainer & options_container) : num_samples(num_samples_in), num_genomic_rate_gc_bias_bins(options_container.getValue<uchar>("number-of-genomic-rate-gc-bias-bins")), noise_rate_priors(num_samples, options_container.getValue<pair<double,double> >("noise-rate-prior")) {

	prng = mt19937(options_container.getValue<int>("random-seed"));
	
	assert(num_samples > 0);

	genomic_count_distributions = genomic_count_distributions_in;

	assert(genomic_count_distributions.size() == num_samples);
	assert(noise_rate_priors.size() == num_samples);

	noise_rates = vector<double>(num_samples, 0);

	for (uint sample_idx = 0; sample_idx < num_samples; sample_idx++) {

		assert(genomic_count_distributions.at(sample_idx).size() == num_genomic_rate_gc_bias_bins);

		noise_rates.at(sample_idx) = sampleGamma(noise_rate_priors.at(sample_idx).first, noise_rate_priors.at(sample_idx).second);
	}

	noise_rate_estimates.push_back(noise_rates);

	genomic_count_log_pmf_cache = vector<vector<vector<vector<double> > > >(num_samples, vector<vector<vector<double> > >(num_genomic_rate_gc_bias_bins, vector<vector<double> >(max_multiplicity + 1, vector<double>(max_kmer_count + 1))));
	noise_count_log_pmf_cache = vector<vector<double> >(num_samples, vector<double>(max_kmer_count + 1));

	updateGenomicCache();
	updateNoiseCache();
}

const vector<vector<NegativeBinomialDistribution> > & CountDistribution::genomicCountDistributions() const {

	return genomic_count_distributions;
}

const vector<double> & CountDistribution::noiseRates() const {

	return noise_rates;
}

const vector<vector<double> > & CountDistribution::noiseRateEstimates() const {

	return noise_rate_estimates;
}

void CountDistribution::setGenomicCountDistributions(const vector<vector<NegativeBinomialDistribution> > & genomic_count_distributions_in) {

	genomic_count_distributions = genomic_count_distributions_in;
}

void CountDistribution::setNoiseRates(const vector<double> & noise_rates_in) {

	noise_rates = noise_rates_in;
}

void CountDistribution::writeNoiseParameterEstimates( const string & output_prefix, const vector<Sample> & samples) const {

	assert(samples.size() == num_samples);
	assert(noise_rates.size() == num_samples);
	assert(noise_rate_priors.size() == num_samples);

    ofstream output_file(output_prefix + "_noise_rate_parameter_estimates.txt");
    assert(output_file.is_open());

    output_file << "#Priors: ";

    output_file << samples.at(0).name << " = (" << noise_rate_priors.at(0).first << ", " << noise_rate_priors.at(0).second << ")";

    for (ushort sample_idx = 1; sample_idx < num_samples; sample_idx++) {

	  	output_file << ", " << samples.at(sample_idx).name << " = (" << noise_rate_priors.at(sample_idx).first << ", " << noise_rate_priors.at(sample_idx).second << ")";
    }

    output_file << "\n#Final estimates: ";

    output_file << samples.at(0).name << " = " << noise_rates.at(0);

    for (ushort sample_idx = 1; sample_idx < num_samples; sample_idx++) {

	  	output_file << ", " << samples.at(sample_idx).name << " = " << noise_rates.at(sample_idx);
    }

    output_file << "\nIteration";

    for (ushort sample_idx = 0; sample_idx < num_samples; sample_idx++) {

	  	output_file << "\t" << samples.at(sample_idx).name;
    }

    output_file << endl;

    for (ushort i = 0; i < noise_rate_estimates.size(); i++) {

		assert(noise_rate_estimates.at(i).size() == num_samples);

    	output_file << i;

    	for (auto & noise_rate_estimate: noise_rate_estimates.at(i)) {

    		output_file << "\t" << noise_rate_estimate;
    	}

	    output_file << endl;
    }

	output_file.close();
}

void CountDistribution::sampleNoiseParameters(const CountAllocation & count_allocation) {

	assert(count_allocation.getCounts().size() == num_samples);

	for (uint sample_idx = 0; sample_idx < num_samples; sample_idx++) {

		assert(count_allocation.getCounts().at(sample_idx).size() == (Utils::uchar_overflow + 1));
		auto count_suff_stats = calcCountSuffStats(count_allocation.getCounts().at(sample_idx));

		noise_rates.at(sample_idx) = sampleGamma(noise_rate_priors.at(sample_idx).first + count_suff_stats.second, noise_rate_priors.at(sample_idx).second/(count_suff_stats.first * noise_rate_priors.at(sample_idx).second + 1));
	}

	noise_rate_estimates.push_back(noise_rates);
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

	assert(genomic_count_log_pmf_cache.size() == num_samples);

	for (ushort sample_idx = 0; sample_idx < num_samples; sample_idx++) {

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

	assert(noise_count_log_pmf_cache.size() == num_samples);

	for (ushort sample_idx = 0; sample_idx < num_samples; sample_idx++) {

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


