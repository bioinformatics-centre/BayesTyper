
/*
HaplotypeFrequencyDistribution.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#include "HaplotypeFrequencyDistribution.hpp"
#include "Utils.hpp"
#include "FrequencyDistribution.hpp"


HaplotypeFrequencyDistribution::HaplotypeFrequencyDistribution() {

	num_haplotype_count = 0;
	num_missing_count = 0;
}

uint HaplotypeFrequencyDistribution::numHaplotypeCount() {

	return num_haplotype_count; 
}

uint HaplotypeFrequencyDistribution::numMissingCount() {

	return num_missing_count; 
}


UniformHaplotypeFrequencyDistribution::UniformHaplotypeFrequencyDistribution(const ushort num_haplotypes) : HaplotypeFrequencyDistribution(), frequency(1 / static_cast<double>(num_haplotypes)) {}

pair<bool, double> UniformHaplotypeFrequencyDistribution::getFrequency(const ushort haplotype_idx) {

	assert(haplotype_idx < Utils::ushort_overflow);
	return make_pair(true, frequency);
}

void UniformHaplotypeFrequencyDistribution::incrementCount(const ushort haplotype_idx) {

	if (haplotype_idx == Utils::ushort_overflow) {

		num_missing_count++;

	} else {

		num_haplotype_count++;
	}
}

void UniformHaplotypeFrequencyDistribution::sampleFrequencies() {

	num_haplotype_count = 0;
	num_missing_count = 0;
}


SparseHaplotypeFrequencyDistribution::SparseHaplotypeFrequencyDistribution(const vector<uint> & non_zero_haplotypes, const uint num_haplotypes, const uint prng_seed) : HaplotypeFrequencyDistribution() {
	
	if (non_zero_haplotypes.empty()) {

		frequency_distribution = new FrequencyDistribution(num_haplotypes, prng_seed);
	
	} else {

		frequency_distribution = new SparseFrequencyDistribution(non_zero_haplotypes.size() / static_cast<double>(num_haplotypes), num_haplotypes, prng_seed);
	}
}

SparseHaplotypeFrequencyDistribution::~SparseHaplotypeFrequencyDistribution() {

	delete frequency_distribution;
}

void SparseHaplotypeFrequencyDistribution::reset() {

	assert(num_haplotype_count == 0);
	assert(num_missing_count == 0);

	frequency_distribution->reset();
}

void SparseHaplotypeFrequencyDistribution::initialize(const vector<uint> & non_zero_haplotypes, const uint num_haplotypes) {

	frequency_distribution->initialize(non_zero_haplotypes);
	frequency_distribution->setSparsity(non_zero_haplotypes.size() / static_cast<double>(num_haplotypes));
}

pair<bool, double> SparseHaplotypeFrequencyDistribution::getFrequency(const ushort haplotype_idx) {

	assert(haplotype_idx < Utils::ushort_overflow);
	return frequency_distribution->getElementFrequency(haplotype_idx);
}

void SparseHaplotypeFrequencyDistribution::incrementCount(const ushort haplotype_idx) {

	if (haplotype_idx == Utils::ushort_overflow) {

		num_missing_count++;

	} else {

		num_haplotype_count++;
		frequency_distribution->incrementObservationCount(haplotype_idx);
	}
}

void SparseHaplotypeFrequencyDistribution::sampleFrequencies() {

	if (num_haplotype_count > 0) {
		
		frequency_distribution->sampleFrequencies(num_haplotype_count);
	}

	num_haplotype_count = 0;
	num_missing_count = 0;
}

