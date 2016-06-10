
/*
DiscreteSampler.cpp - This file is part of BayesTyper (v0.9)


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


#include <limits>

#include "Utils.hpp"
#include "DiscreteSampler.hpp"

DiscreteSampler::DiscreteSampler() {}

DiscreteSampler::DiscreteSampler(uint size_guess) {

	has_finite = false;
	cum_probs.reserve(size_guess);
}

DiscreteSampler::DiscreteSampler(vector<double> & cum_probs_in) {

	cum_probs = cum_probs_in;
}

void DiscreteSampler::addOutcome(double prob) {

	assert(std::isfinite(prob));

	if (cum_probs.size() == 1) {

		cum_probs.push_back(prob);
	
	} else {

		cum_probs.push_back(prob + cum_probs.back());
	}
}

uint DiscreteSampler::search(double random01_sample_scaled) {

	if (cum_probs.size() > 1) {

		auto sample_iter = upper_bound(cum_probs.begin(), cum_probs.end(), random01_sample_scaled);
		uint sample_idx = sample_iter - cum_probs.begin();
		assert(sample_iter != cum_probs.end());
		assert(sample_idx >= 0);

		return sample_idx;

	} else {

		assert(cum_probs.size() == 1);
		return 0;
	}
}

uint DiscreteSampler::sample(mt19937 & prng) {

	double random01_sample_scaled = generate_canonical<double,std::numeric_limits<double>::digits>(prng) * cum_probs.back();

	return search(random01_sample_scaled);
}

double DiscreteSampler::probSum() {

	assert(!(cum_probs.empty()));
	return cum_probs.back();
}

void DiscreteSampler::clear() {

	cum_probs.clear();
}

LogDiscreteSampler::LogDiscreteSampler() {}

LogDiscreteSampler::LogDiscreteSampler(uint size_guess) : DiscreteSampler(size_guess) {}

LogDiscreteSampler::LogDiscreteSampler(vector<double> & cum_log_probs_in) : DiscreteSampler(cum_log_probs_in) {}

void LogDiscreteSampler::addOutcome(double log_prob) {

	if (std::isfinite(log_prob)) {

		has_finite = true;

		if (cum_probs.size() == 0) {

			cum_probs.push_back(log_prob);

		} else {

			double log_sum = Utils::logAddition(log_prob, cum_probs.back());
			assert(std::isfinite(log_sum));
			cum_probs.push_back(log_sum); 
		}

	} else {

		assert(log_prob == -numeric_limits<double>::infinity());

		if (cum_probs.size() == 0) {

			cum_probs.push_back(numeric_limits<double>::lowest());
		
		} else {
			
			cum_probs.push_back(cum_probs.back());
		}				
	}
}

uint LogDiscreteSampler::sample(mt19937 & prng) {

	if (has_finite) {

		double random01_log_sample_scaled = log(generate_canonical<double,std::numeric_limits<double>::digits>(prng)) + cum_probs.back();
		return search(random01_log_sample_scaled);

	} else {

		std::uniform_int_distribution<int> int_dist(0, cum_probs.size() - 1);
		return int_dist(prng);
	}
}

