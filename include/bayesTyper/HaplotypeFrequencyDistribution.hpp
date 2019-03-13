
/*
HaplotypeFrequencyDistribution.hpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#ifndef __bayesTyper__HaplotypeFrequencyDistribution_hpp
#define __bayesTyper__HaplotypeFrequencyDistribution_hpp

#include "Utils.hpp"
#include "FrequencyDistribution.hpp"

class HaplotypeFrequencyDistribution {

	public:

		HaplotypeFrequencyDistribution();
		virtual ~HaplotypeFrequencyDistribution() {};

		uint numHaplotypeCount();
		uint numMissingCount();

		virtual void reset() = 0;
		virtual void initialize(const vector<uint> &, const uint) = 0;

		virtual pair<bool, double> getFrequency(const ushort) = 0;
		virtual void incrementCount(const ushort) = 0;

		virtual void sampleFrequencies() = 0;

	protected:

		uint num_haplotype_count;
		uint num_missing_count;
};

class UniformHaplotypeFrequencyDistribution : public HaplotypeFrequencyDistribution {

	public:

		UniformHaplotypeFrequencyDistribution(const ushort);
		~UniformHaplotypeFrequencyDistribution() {};

		void reset() {};
		void initialize(const vector<uint> &, const uint) {};

		pair<bool, double> getFrequency(const ushort);
		void incrementCount(const ushort);

		void sampleFrequencies();

	private:

		const double frequency;
};

class SparseHaplotypeFrequencyDistribution : public HaplotypeFrequencyDistribution {

	public:

		SparseHaplotypeFrequencyDistribution(const vector<uint> &, const uint, const uint);
		~SparseHaplotypeFrequencyDistribution();

		void reset();
		void initialize(const vector<uint> &, const uint);

		pair<bool, double> getFrequency(const ushort);
		void incrementCount(const ushort);

		void sampleFrequencies();

	protected:

		FrequencyDistribution * frequency_distribution;
};

#endif	

