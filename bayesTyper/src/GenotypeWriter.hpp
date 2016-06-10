
/*
GenotypeWriter.hpp - This file is part of BayesTyper (v0.9)


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


using namespace std;

class GenotypeWriter {
    
	private:

		const ushort num_samples;

		const string output_prefix;
		const double num_gibbs_samples;
		const ushort min_observed_kmers;

		const ChromosomePloidy chromosome_ploidy;

    	ProducerConsumerQueue<vector<Genotypes *> * > * genotype_queue;
		vector<thread> writer_threads; 

		struct AlleleCounts {

			vector<uint> counts;
			uint num_alleles;

			AlleleCounts(const uint num_alternative_alleles) : counts(num_alternative_alleles, 0), num_alleles(0) {}

			void addAllele(const ushort allele_idx) {

			    num_alleles++;

			    if (allele_idx > 0) {
			        
			        counts.at(allele_idx - 1)++;
			    }
			}
		};

		void writeTemporaryFile();

		uint genotypeToOneDimensionalIndex(const pair<ushort, ushort>, const ushort);

		vector<vector<string> > getAlleleFilterFlags(const VariantKmerStats &, const VariantInfo &);
		void estimateGenotypesAndAlleleCounts(const PosteriorContainer &, const VariantInfo &, vector<vector<string> > *, vector<vector<float> > *, vector<vector<float> > *, AlleleCounts *, vector<vector<string> > *); 
		vector<float> calculateAlleleCallProbabilities(const vector<vector<float> > &, const VariantInfo &);
		
		void writeQualityAndFilter(ofstream *, Utils::FilterStatus, const vector<float> &, const VariantInfo &);
		void writeAlleleCounts(ofstream *, const AlleleCounts &);
		void writeAlleleCallProbabilities(ofstream *, const vector<float> &);
		void writeSamples(ofstream *, const vector<vector<string> > &, const vector<vector<float> > &, const vector<vector<float> > &, const VariantKmerStats &, const VariantInfo &, vector<vector<string> > *);

		template <typename ValueType>
		void writeAlleleFormatField(ofstream *, const vector<ValueType> &, const bool);
	
		template <typename FileType>
		uint parseAndWriteGenotypes(FileType *,  const Regions &, const vector<Sample> &, const string &, const string &, const uint);

		string writeHeader(const vector<Sample> &, const string &, const string &, const string &);		
		void writeGenotypesVariant(ofstream *, string *, const string &);
		
		void writeUnsupportedVariant(ofstream *, const string &, const string &, const string &, const uint);

	public:

		GenotypeWriter(const vector<Sample> &, const OptionsContainer &);

		void addGenotypes(vector<Genotypes *> *);
		void completedGenotyping();
		uint writeGenotypesToVariantFile(const string &, const Regions &, const vector<Sample> &, const string, const string &, const uint);	
};

#endif