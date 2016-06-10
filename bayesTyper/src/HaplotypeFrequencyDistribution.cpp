
/*
HaplotypeFrequencyDistribution.cpp - This file is part of BayesTyper (v0.9)


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


HaplotypeFrequencyDistribution::HaplotypeFrequencyDistribution(FrequencyDistribution * frequency_distribution_in) : frequency_distribution(frequency_distribution_in) {

    sum_haplotype_count = 0;
    sum_missing_count = 0;
}


HaplotypeFrequencyDistribution::~HaplotypeFrequencyDistribution() {

	delete frequency_distribution;
}


void HaplotypeFrequencyDistribution::reset() {

	assert(sum_haplotype_count == 0);
	assert(sum_missing_count == 0);

	frequency_distribution->reset();
}

pair<bool, double> HaplotypeFrequencyDistribution::getElementFrequency(const ushort element_idx) {

	assert(element_idx < Utils::ushort_overflow);
	return frequency_distribution->getElementFrequency(element_idx);
}


void HaplotypeFrequencyDistribution::incrementObservationCount(const ushort element_idx) {

	if (element_idx == Utils::ushort_overflow) {

		sum_missing_count++;

	} else {

		sum_haplotype_count++;
		frequency_distribution->incrementObservationCount(element_idx);
	}
}


uint HaplotypeFrequencyDistribution::sumHaplotypeCount() {

	return sum_haplotype_count; 
}


uint HaplotypeFrequencyDistribution::sumMissingCount() {

	return sum_missing_count; 
}


void HaplotypeFrequencyDistribution::sampleFrequencies() {

	if (sum_haplotype_count > 0) {
		
		frequency_distribution->sampleFrequencies(sum_haplotype_count);
	}

	sum_haplotype_count = 0;
	sum_missing_count = 0;
}

