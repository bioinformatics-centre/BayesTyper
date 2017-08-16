
/*
GenotypeWriter.hpp - This file is part of BayesTyper (v1.1)


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


#ifndef __bayesTyper__GenotypeWriter_hpp
#define __bayesTyper__GenotypeWriter_hpp

#include <string>
#include <map>
#include <vector>
#include <unordered_map>
#include <thread>
#include <iostream>
#include <fstream>

#include "boost/functional/hash.hpp"
#include "ProducerConsumerQueue.hpp"

#include "Utils.hpp"
#include "Genotypes.hpp"
#include "Sample.hpp"
#include "OptionsContainer.hpp"
#include "ChromosomePloidy.hpp"
#include "Regions.hpp"
#include "KmerStats.hpp"


using namespace std;

class GenotypeWriter {
    
	private:

		const vector<Sample> & samples;
		const string output_prefix;

    	ProducerConsumerQueue<vector<Genotypes *> * > * genotypes_queue;
		vector<thread> writer_threads; 

		vector<vector<ulong> > sample_fak_estimates;
		vector<vector<ulong> > sample_mac_estimates;

		void writeTemporaryFile();

		template <typename ValueType>
		void writeAlleleField(ofstream *, const vector<ValueType> &);

		void writeQualityAndFilter(ofstream *, const Genotypes::VariantStats &, const VariantInfo &);
		void writeVariantStats(ofstream *, const Genotypes::VariantStats &, const VariantInfo &);
		void writeAlleleCover(ofstream *, vector<ushort> *, VariantInfo *);
		void writeSamples(ofstream *, const vector<Genotypes::SampleStats> &, const VariantInfo &);
		void writeAlleleKmerStats(ofstream *, const AlleleKmerStats &);

		void addAlleleKmerEstimates(const vector<Genotypes::SampleStats> &);

		void writeSampleAttributeCumDistFunc(const string &, const vector<vector<ulong> > &, const pair<uint, uint> &);
	
		template <typename FileType>
		uint parseAndWriteGenotypes(FileType *,  const Regions &, const string &, const string &, const uint);

		string writeHeader(const vector<Sample> &, const string &, const string &);		
		void writeUnsupportedVariant(ofstream *, const string &, const string &, const uint, const ChromosomePloidy &);

	public:

		GenotypeWriter(const vector<Sample> &, const string &, const ushort);

		void addGenotypes(vector<Genotypes *> *);
		void completedGenotyping();

		void writeSampleAlleleKmerFractionCumDistFunc();
		void writeSampleAlleleKmerCoverageCumDistFunc();

		uint writeGenotypesToVariantCallFormat(const string &, const Regions &, const string &, const string &, const uint);	
};

#endif