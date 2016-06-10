
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

static const ushort num_noise_parameter_gibbs_samples = 10000;

CountDistribution::CountDistribution(const ushort num_samples_in, const ushort num_noise_sources_in, const vector<NegativeBinomialDistribution> & genomic_count_distributions_in, const OptionsContainer & options_container) : num_samples(num_samples_in), num_noise_sources(num_noise_sources_in), noise_zero_inflation_priors(num_samples, options_container.getValue<vector<pair<double,double> > >("noise-zero-inflation-priors")) , noise_rate_priors(num_samples, options_container.getValue<vector<pair<double,double> > >("noise-rate-priors")) {

	prng = mt19937(options_container.getValue<int>("random-seed"));
	
	assert(num_samples > 0);
	assert(num_noise_sources > 0);
	assert(num_noise_sources <= 8);

	genomic_count_distributions = genomic_count_distributions_in;

	assert(genomic_count_distributions.size() == num_samples);
	assert(noise_zero_inflation_priors.size() == num_samples);
	assert(noise_rate_priors.size() == num_samples);

	for (uint sample_idx = 0; sample_idx < num_samples; sample_idx++) {

		assert(noise_zero_inflation_priors.at(sample_idx).size() == num_noise_sources);
		assert(noise_rate_priors.at(sample_idx).size() == num_noise_sources);
	}

	noise_zero_inflations = vector<vector<double> >(num_samples, vector<double>(num_noise_sources));
	noise_rates = vector<vector<double> >(num_samples, vector<double>(num_noise_sources));

	for (uint sample_idx = 0; sample_idx < num_samples; sample_idx++) {

		for (uint noise_idx = 0; noise_idx < num_noise_sources; noise_idx++) {

			noise_zero_inflations.at(sample_idx).at(noise_idx) = sampleBeta(noise_zero_inflation_priors.at(sample_idx).at(noise_idx).first, noise_zero_inflation_priors.at(sample_idx).at(noise_idx).second);
			noise_rates.at(sample_idx).at(noise_idx) = sampleGamma(noise_rate_priors.at(sample_idx).at(noise_idx).first, noise_rate_priors.at(sample_idx).at(noise_idx).second);
		}
	}

	max_multiplicity_cache = 2;
	max_kmer_count = Utils::uchar_overflow;

	genomic_count_log_pmf_cache = vector<vector<vector<double> > >(num_samples, vector<vector<double> >(max_multiplicity_cache + 1, vector<double>(max_kmer_count + 1)));
	combined_noise_count_log_prob_cache = vector<vector<double> >(num_samples, vector<double>(max_kmer_count + 1));
	combined_genotype_noise_count_log_prob_cache = vector<vector<vector<double> > >(num_samples, vector<vector<double> >(max_multiplicity_cache + 1, vector<double>(max_kmer_count + 1)));

	updateCache();
}

const vector<NegativeBinomialDistribution> & CountDistribution::genomicCountDistributions() const {

	return genomic_count_distributions;
}

const vector<vector<double> > & CountDistribution::noiseZeroInflations() const {

	return noise_zero_inflations;
}

const vector<vector<double> > & CountDistribution::noiseRates() const {

	return noise_rates;
}

const vector<vector<vector<double> > > & CountDistribution::noiseZeroInflationsSamples() const {

	return noise_zero_inflation_samples;
}

const vector<vector<vector<double> > > & CountDistribution::noiseRatesSamples() const {

	return noise_rate_samples;
}

void CountDistribution::setGenomicCountDistributions(const vector<NegativeBinomialDistribution> & genomic_count_distributions_in) {

	genomic_count_distributions = genomic_count_distributions_in;
}

void CountDistribution::setNoiseZeroInflations(const vector<vector<double> > & noise_zero_inflations_in) {

	noise_zero_inflations = noise_zero_inflations_in;
}

void CountDistribution::setNoiseRates(const vector<vector<double> > & noise_rates_in) {

	noise_rates = noise_rates_in;
}

void CountDistribution::updateMaxMultiplicityCache(const uchar max_multiplicity_cache_in) {

	assert(max_multiplicity_cache_in <= Utils::uchar_overflow);
	max_multiplicity_cache = max_multiplicity_cache_in;

	genomic_count_log_pmf_cache = vector<vector<vector<double> > >(num_samples, vector<vector<double> >(max_multiplicity_cache + 1, vector<double>(max_kmer_count + 1)));
	combined_genotype_noise_count_log_prob_cache = vector<vector<vector<double> > >(num_samples, vector<vector<double> >(max_multiplicity_cache + 1, vector<double>(max_kmer_count + 1)));

	updateCache();
}

void CountDistribution::writeNoiseSamples(const vector<vector<vector<double> > > & noise_samples, const string & filename, const ushort noise_idx, const vector<vector<pair<double,double> > > & noise_priors, const vector<vector<double> > & final_estimates, const vector<Sample> & samples) const {

	assert(final_estimates.size() == num_samples);
	assert(noise_priors.size() == num_samples);
	assert(samples.size() == num_samples);

    ofstream output_file(filename);
    assert(output_file.is_open());

    output_file << "#Priors: ";

    output_file << samples.at(0).name << " = (" << noise_priors.at(0).at(noise_idx).first << ", " << noise_priors.at(0).at(noise_idx).second << ")";

    for (ushort sample_idx = 1; sample_idx < num_samples; sample_idx++) {

	  	output_file << ", " << samples.at(sample_idx).name << " = (" << noise_priors.at(sample_idx).at(noise_idx).first << ", " << noise_priors.at(sample_idx).at(noise_idx).second << ")";
    }

    output_file << "\n#Final estimates: ";

    output_file << samples.at(0).name << " = " << final_estimates.at(0).at(noise_idx);

    for (ushort sample_idx = 1; sample_idx < num_samples; sample_idx++) {

	  	output_file << ", " << samples.at(sample_idx).name << " = " << final_estimates.at(sample_idx).at(noise_idx);
    }

    output_file << "\nIteration";

    for (ushort sample_idx = 0; sample_idx < num_samples; sample_idx++) {

	  	output_file << "\t" << samples.at(sample_idx).name;
    }

    output_file << endl;

    ushort iteration = 1;

    for (auto & noise_sample: noise_samples) {

		assert(noise_sample.size() == num_samples);

    	output_file << iteration;

    	for (auto & value: noise_sample) {

    		output_file << "\t" << value.at(noise_idx);
    	}

	    output_file << endl;
    	iteration += 1;
    }

	output_file.close();
}

void CountDistribution::writeParameterSamples(const vector<Sample> & samples) const {

	for (ushort i = 0; i < num_noise_sources; i++) {

    	writeNoiseSamples(noise_zero_inflation_samples, "model_noise" + to_string(i + 1) + "_zero_inflation_parameter_estimations.txt", i, noise_zero_inflation_priors, noise_zero_inflations, samples);
    	writeNoiseSamples(noise_rate_samples, "model_noise" + to_string(i + 1) + "_rate_parameter_estimations.txt", i, noise_rate_priors, noise_rates, samples);
	}
}

double CountDistribution::sampleGamma(const double shape, const double scale) {

	assert(shape > 0);
	assert(scale > 0);

	gamma_dist.param(gamma_distribution<>::param_type(shape, scale));
	double sample = gamma_dist(prng);

	assert(sample > 0);

	return sample;
}

double CountDistribution::sampleBeta(const double alpha, const double beta) {

	assert(alpha > 0);
	assert(beta > 0);

	double gamma_sample_alpha = sampleGamma(alpha, 1);
	double gamma_sample_beta = sampleGamma(beta, 1);

	double sample = gamma_sample_alpha/(gamma_sample_alpha + gamma_sample_beta);

	assert(sample > 0);
	assert(sample < 1);

	return sample;
}

void CountDistribution::sampleParameters(const CountAllocation & count_allocation) {

	sampleNoiseParameters(count_allocation.noiseCounts());

	noise_zero_inflation_samples.push_back(noise_zero_inflations);
	noise_rate_samples.push_back(noise_rates);

	updateCache();
}

void CountDistribution::sampleNoiseParameters(const vector<vector<unordered_map<uchar,ulong> > > & noise_counts) {

	for (uint sample_idx = 0; sample_idx < num_samples; sample_idx++) {

		for (uint noise_idx = 0; noise_idx < num_noise_sources; noise_idx++) {
				
			auto count_suff_stats = calcCountSuffStats(noise_counts.at(sample_idx).at(noise_idx));
			auto noise_zero_count_find_result = noise_counts.at(sample_idx).at(noise_idx).find(0);

			for (uint i = 0; i < num_noise_parameter_gibbs_samples; i++) {

				ulong num_poisson_observations = count_suff_stats.first;

				if (noise_zero_count_find_result != noise_counts.at(sample_idx).at(noise_idx).end()) {

					double poisson_zero_prob = (1 - noise_zero_inflations.at(sample_idx).at(noise_idx))*exp(-noise_rates.at(sample_idx).at(noise_idx));
					double poisson_indicator_posterior = poisson_zero_prob/(noise_zero_inflations.at(sample_idx).at(noise_idx) + poisson_zero_prob);

					assert(poisson_indicator_posterior > 0);
					assert(poisson_indicator_posterior < 1);

					binomial_dist.param(binomial_distribution<>::param_type(noise_zero_count_find_result->second, poisson_indicator_posterior));
					num_poisson_observations = count_suff_stats.first - noise_zero_count_find_result->second + binomial_dist(prng);
				}

				assert(num_poisson_observations <= count_suff_stats.first);

				noise_rates.at(sample_idx).at(noise_idx) = sampleGamma(noise_rate_priors.at(sample_idx).at(noise_idx).first + count_suff_stats.second, noise_rate_priors.at(sample_idx).at(noise_idx).second/(num_poisson_observations*noise_rate_priors.at(sample_idx).at(noise_idx).second + 1));
				noise_zero_inflations.at(sample_idx).at(noise_idx) = sampleBeta(noise_zero_inflation_priors.at(sample_idx).at(noise_idx).first + count_suff_stats.first - num_poisson_observations, noise_zero_inflation_priors.at(sample_idx).at(noise_idx).second + num_poisson_observations);
			}
		}
	}
}

double CountDistribution::calcCombinedGenotypeNoiseCountLogProb(const ushort sample_idx, const uchar genotype_kmer_multiplicity, const uchar kmer_count) const {

	return combined_genotype_noise_count_log_prob_cache.at(sample_idx).at(genotype_kmer_multiplicity).at(kmer_count);		
}

double CountDistribution::calcCombinedGenotypeNoiseCountLogProbInternal(const ushort sample_idx, const uchar genotype_kmer_multiplicity, const uchar kmer_count) const {

	auto genotype_count_log_posterior = calcGenotypeCountLogPosterior(sample_idx, genotype_kmer_multiplicity, kmer_count);
	return genotype_count_log_posterior->probSum();
}

double CountDistribution::calcCombinedNoiseCountLogProb(const ushort sample_idx, const uchar kmer_count) const {

	return combined_noise_count_log_prob_cache.at(sample_idx).at(kmer_count);
}

double CountDistribution::calcCombinedNoiseCountLogProbInternal(const ushort sample_idx, const vector<uchar> & noise_indices, const uchar kmer_count) const {

	assert(!(noise_indices.empty()));
	assert(noise_indices.size() <= num_noise_sources);

	double kmer_count_log_prob = -std::numeric_limits<double>::infinity();
	auto zero_inflation_combinations = Combinator::enumerateBinaryCombinations(noise_indices.size());

	for (auto & zero_inflation_combination: zero_inflation_combinations) {

		assert(zero_inflation_combination.size() == noise_indices.size());

		double noise_zero_inflation_log_prob = 0;
		double noise_convoluted_rate = 0;

		bool has_noise_rate = false;

		for (ushort i = 0; i < zero_inflation_combination.size(); i++) {

			if (zero_inflation_combination.at(i)) {

				noise_zero_inflation_log_prob += log(noise_zero_inflations.at(sample_idx).at(noise_indices.at(i)));

			} else {

				has_noise_rate = true;

				noise_zero_inflation_log_prob += log(1 - noise_zero_inflations.at(sample_idx).at(noise_indices.at(i)));
				noise_convoluted_rate += noise_rates.at(sample_idx).at(noise_indices.at(i));
			}
		}

		if (has_noise_rate) {

			assert(noise_convoluted_rate > Utils::double_underflow);
			kmer_count_log_prob = Utils::logAddition(noise_zero_inflation_log_prob + poissonLogPmf(noise_convoluted_rate, kmer_count), kmer_count_log_prob);

		} else if (kmer_count == 0) {

			assert(Utils::doubleCompare(noise_convoluted_rate, 0));
			kmer_count_log_prob = Utils::logAddition(noise_zero_inflation_log_prob, kmer_count_log_prob);
		}
	}

	assert(kmer_count_log_prob <= 0);

	return kmer_count_log_prob;
}

unique_ptr<LogDiscreteSampler> CountDistribution::calcGenotypeCountLogPosterior(const ushort sample_idx, const uchar genotype_kmer_multiplicity, const uchar kmer_count) const {

	unique_ptr<LogDiscreteSampler> genotype_count_log_posterior(new LogDiscreteSampler(kmer_count + 1));

	for (uint genomic_count = 0; genomic_count <= kmer_count; genomic_count++) {

		double genotype_count_log_prob = genomicCountLogPmf(sample_idx, genotype_kmer_multiplicity, genomic_count) + calcCombinedNoiseCountLogProb(sample_idx, kmer_count - genomic_count);
		genotype_count_log_posterior->addOutcome(genotype_count_log_prob);
	}

	return genotype_count_log_posterior;
}

unique_ptr<LogDiscreteSampler> CountDistribution::calcNoiseSplitCountLogPosterior(const ushort sample_idx, const uchar split_noise_idx, const uchar kmer_count) const {

	unique_ptr<LogDiscreteSampler> noise_split_count_log_posterior(new LogDiscreteSampler(kmer_count + 1));

	assert(split_noise_idx < (num_noise_sources - 1));

	vector<uchar> noise_indices(num_noise_sources - split_noise_idx - 1); 
	iota(noise_indices.begin(), noise_indices.end(), split_noise_idx + 1);

	assert(!(noise_indices.empty()));

	for (uint split_count = 0; split_count <= kmer_count; split_count++) {

		double noise_split_count_log_prob = noiseCountLogPmf(sample_idx, split_noise_idx, split_count) + calcCombinedNoiseCountLogProbInternal(sample_idx, noise_indices, kmer_count - split_count);
		noise_split_count_log_posterior->addOutcome(noise_split_count_log_prob);
	}

	return noise_split_count_log_posterior;
}

void CountDistribution::updateCache() {

	for (uint i = 0; i < num_samples; i++) {

		for (uint k = 0; k <= max_multiplicity_cache; k++) {

			for (uint kmer_count = 0; kmer_count <= max_kmer_count; kmer_count++) {

				genomic_count_log_pmf_cache.at(i).at(k).at(kmer_count) = genomicCountLogPmfInternal(i, k, kmer_count);
			}
		}
	}

	vector<uchar> noise_indices(num_noise_sources); 
	iota(noise_indices.begin(), noise_indices.end(), 0);	

	for (uint i = 0; i < num_samples; i++) {

		for (uint kmer_count = 0; kmer_count <= max_kmer_count; kmer_count++) {

			combined_noise_count_log_prob_cache.at(i).at(kmer_count) = calcCombinedNoiseCountLogProbInternal(i, noise_indices, kmer_count);
		}
	}

	for (uint i = 0; i < num_samples; i++) {

		for (uint k = 0; k <= max_multiplicity_cache; k++) {

			for (uint kmer_count = 0; kmer_count <= max_kmer_count; kmer_count++) {

				combined_genotype_noise_count_log_prob_cache.at(i).at(k).at(kmer_count) = calcCombinedGenotypeNoiseCountLogProbInternal(i, k, kmer_count);
			}
		}
	}
}

double CountDistribution::genomicCountLogPmf(const ushort sample_idx, const uchar genotype_kmer_multiplicity, const uchar kmer_count) const {

	return genomic_count_log_pmf_cache.at(sample_idx).at(genotype_kmer_multiplicity).at(kmer_count);
}

double CountDistribution::genomicCountLogPmfInternal(const ushort sample_idx, const uchar genotype_kmer_multiplicity, const uchar kmer_count) const {

	if (genotype_kmer_multiplicity == 0) {

		if (kmer_count == 0) {

			return 0;

		} else {

			return -numeric_limits<double>::infinity();
		}

	} else if (kmer_count == max_kmer_count) {

		LogDiscreteSampler lower_kmer_count_log_prob = LogDiscreteSampler(kmer_count);

		for (uint lower_kmer_count = 0; lower_kmer_count < kmer_count; lower_kmer_count++) {

			lower_kmer_count_log_prob.addOutcome(genomic_count_log_pmf_cache.at(sample_idx).at(genotype_kmer_multiplicity).at(lower_kmer_count));
		}

		double lower_kmer_count_sum_log = lower_kmer_count_log_prob.probSum();

		if (lower_kmer_count_sum_log < 0) {

			return log(1 - exp(lower_kmer_count_sum_log));

		} else {

			return -numeric_limits<double>::infinity();
		}

	} else {

		return genomic_count_distributions.at(sample_idx).logPmf(kmer_count, genotype_kmer_multiplicity);
	}

}

double CountDistribution::poissonLogPmf(const double rate, const uchar count) const {

	assert(rate > Utils::double_underflow);
	return count * log(rate) - rate - boost::math::lgamma(count + 1);
}

double CountDistribution::noiseCountLogPmf(const ushort sample_idx, const uchar noise_idx, const uchar kmer_count) const {

	if (kmer_count == 0) {

		return Utils::logAddition(log(noise_zero_inflations.at(sample_idx).at(noise_idx)), log(1-noise_zero_inflations.at(sample_idx).at(noise_idx)) + poissonLogPmf(noise_rates.at(sample_idx).at(noise_idx), kmer_count));

	} else {

		return log(1-noise_zero_inflations.at(sample_idx).at(noise_idx)) + poissonLogPmf(noise_rates.at(sample_idx).at(noise_idx), kmer_count);
	}
}

pair<ulong,ulong> CountDistribution::calcCountSuffStats(const unordered_map<uchar,ulong> & counts) const {

	ulong num_observations = 0;
	ulong count_sum = 0;

	for (auto & count : counts) {

		num_observations += count.second;
		count_sum += count.first * count.second;
	}

	return pair<ulong,ulong> (num_observations, count_sum);
}
