
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


static const uint max_nb_kmer_multiplicity = 32;
static const uint min_nb_kmer_count = 10000;

static const uchar num_genomic_rate_gc_bias_bins = 1;
static const uchar max_kmer_multiplicity = Utils::uchar_overflow;
static const uchar max_kmer_count = Utils::uchar_overflow;


CountDistribution::CountDistribution(const vector<Sample> & samples_in, const OptionsContainer & options_container) : samples(samples_in), noise_rate_priors(samples.size(), options_container.getValue<pair<float,float> >("noise-rate-prior")) {

	prng = mt19937(options_container.getValue<uint>("random-seed"));
	
	genomic_count_distributions = vector<vector<NegativeBinomialDistribution> >(samples.size(), vector<NegativeBinomialDistribution>(num_genomic_rate_gc_bias_bins));

	noise_rates = vector<double>(samples.size(), 0);
	genomic_count_log_pmf_cache = vector<vector<vector<vector<double> > > >(samples.size(), vector<vector<vector<double> > >(num_genomic_rate_gc_bias_bins, vector<vector<double> >(max_kmer_multiplicity + 1, vector<double>(max_kmer_count + 1))));
	noise_count_log_pmf_cache = vector<vector<double> >(samples.size(), vector<double>(max_kmer_count + 1));

	resetNoiseRates();
	updateGenomicCache();
}


void CountDistribution::setGenomicCountDistributions(const vector<vector<vector<KmerStats> > > & intercluster_kmer_stats, const string & output_prefix) {

    cout << "[" << Utils::getLocalTime() << "] Estimating genomic haploid kmer count distribution(s) from parameter kmers ...\n" << endl;

    ofstream genomic_outfile(output_prefix + ".txt");

    if (!genomic_outfile.is_open()) {

        cerr << "\nERROR: Unable to write file " << output_prefix + ".txt" << "\n" << endl;
        exit(1);
    }

    genomic_outfile << "Sample\tMean\tVariance" << endl;

	assert(intercluster_kmer_stats.size() == samples.size());
	assert(intercluster_kmer_stats.size() == genomic_count_distributions.size());

    assert(num_genomic_rate_gc_bias_bins == 1);
    assert(max_nb_kmer_multiplicity <= Utils::uchar_overflow);

    for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

		assert(intercluster_kmer_stats.at(sample_idx).size() == num_genomic_rate_gc_bias_bins);    	
		assert(intercluster_kmer_stats.at(sample_idx).size() == genomic_count_distributions.at(sample_idx).size());

        for (ushort bias_idx = 0; bias_idx < num_genomic_rate_gc_bias_bins; bias_idx++) {

			assert(intercluster_kmer_stats.at(sample_idx).at(bias_idx).size() == (Utils::uchar_overflow + 1));

        	uint max_genomic_kmer_count = 0;
        	ushort max_genomic_kmer_multiplicity = 0;

			for (ushort kmer_genomic_multiplicity = 1; kmer_genomic_multiplicity <= max_nb_kmer_multiplicity; kmer_genomic_multiplicity++) {

				auto kmer_genomic_count = intercluster_kmer_stats.at(sample_idx).at(bias_idx).at(kmer_genomic_multiplicity).getCount();

				if (kmer_genomic_count > max_genomic_kmer_count) {

					max_genomic_kmer_count = kmer_genomic_count;
					max_genomic_kmer_multiplicity = kmer_genomic_multiplicity;
				}
        	}

			if (max_genomic_kmer_count < min_nb_kmer_count) {

				cout << "\nWARNING: Low number of kmers used for negative binomial parameters estimation for sample " << samples.at(sample_idx).name << " (" << max_genomic_kmer_count << " < " << min_nb_kmer_count << ")" << endl;
				cout << "WARNING: The mean and variance estimates might be biased due to the genome used being too small, too variant dense and/or too repetitive\n" << endl;
			}

			assert(max_genomic_kmer_multiplicity > 0);

    		auto kmer_mean = intercluster_kmer_stats.at(sample_idx).at(bias_idx).at(max_genomic_kmer_multiplicity).getMean();
   	        assert(kmer_mean.second);

    		auto kmer_var = intercluster_kmer_stats.at(sample_idx).at(bias_idx).at(max_genomic_kmer_multiplicity).getVariance();
    		assert(kmer_var.second);

	        auto nb_para = NegativeBinomialDistribution::momentsToParameters(kmer_mean.first, kmer_var.first);
	        nb_para.second /= max_genomic_kmer_multiplicity;

        	genomic_count_distributions.at(sample_idx).at(bias_idx).setParameters(nb_para);

        	auto nb_mean = genomic_count_distributions.at(sample_idx).at(bias_idx).mean();
        	auto nb_var = genomic_count_distributions.at(sample_idx).at(bias_idx).var();

            cout << "[" << Utils::getLocalTime() << "] Estimated negative binomial (mean = " << nb_mean << ", var = " << nb_var << ") for sample " << samples.at(sample_idx).name << " using " << max_genomic_kmer_count << " parameter kmers (multiplicity = " << max_genomic_kmer_multiplicity << ")" << endl;

	        genomic_outfile << samples.at(sample_idx).name << "\t" << nb_mean << "\t" << nb_var << endl;
        }
    }

    genomic_outfile.close();
    updateGenomicCache();

	cout << "\n[" << Utils::getLocalTime() << "] Wrote genomic parameters to " << output_prefix << ".txt" << endl;
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
		count_sum += i * counts.at(i);
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

			assert(genomic_count_log_pmf_cache.at(sample_idx).at(bias_idx).size() == static_cast<uint>(max_kmer_multiplicity + 1));

			for (ushort kmer_multiplicity = 0; kmer_multiplicity <= max_kmer_multiplicity; kmer_multiplicity++) {

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

		} while (!Utils::doubleCompare(prev_genomic_count_log_pmf, genomic_count_log_pmf));
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

		} while (!Utils::doubleCompare(prev_noise_count_log_pmf, noise_count_log_pmf));
	} 

	assert(isfinite(noise_count_log_pmf));
	assert(noise_count_log_pmf <= 0);

	return noise_count_log_pmf;
}

double CountDistribution::poissonLogProb(const uint value, const double rate) const {

	return value * log(rate) - rate - boost::math::lgamma(value + 1);
}


