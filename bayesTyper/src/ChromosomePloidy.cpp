
/*
ChromosomePloidy.cpp - This file is part of BayesTyper (v0.9)


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


#include <vector>

#include "ChromosomePloidy.hpp"
#include "Utils.hpp"
#include "Sample.hpp"

ChromosomePloidy::ChromosomePloidy(const vector<Sample> & samples) {

	autosomal.reserve(samples.size());
	chrX.reserve(samples.size());
	chrY.reserve(samples.size());
	
	for (auto &sample: samples) {

		autosomal.push_back(Utils::Ploidy::Diploid);

		if (sample.sex == Utils::Sex::Male) {

			chrX.push_back(Utils::Ploidy::Haploid);
			chrY.push_back(Utils::Ploidy::Haploid);
		
		} else {

			assert(sample.sex == Utils::Sex::Female);

			chrX.push_back(Utils::Ploidy::Diploid);
			chrY.push_back(Utils::Ploidy::Null);			
		}
	} 
}


const vector<Utils::Ploidy> & ChromosomePloidy::getPloidy(const Utils::ChromosomeClass chromosome_class) const {

	if (chromosome_class == Utils::ChromosomeClass::Autosomal) {

		return autosomal;
	
	} else if (chromosome_class == Utils::ChromosomeClass::X) {

		return chrX;
	
	} else {

		assert(chromosome_class == Utils::ChromosomeClass::Y);

		return chrY;
	}
}