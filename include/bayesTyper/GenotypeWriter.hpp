
/*
GenotypeWriter.hpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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

#include "boost/iostreams/filtering_stream.hpp"
#include "boost/iostreams/filter/gzip.hpp"

#include "Utils.hpp"
#include "Genotypes.hpp"
#include "Sample.hpp"
#include "OptionsContainer.hpp"
#include "ChromosomePloidy.hpp"
#include "Regions.hpp"
#include "KmerStats.hpp"
#include "Chromosomes.hpp"
#include "ProducerConsumerQueue.hpp"
#include "Filters.hpp"


using namespace std;

class GenotypeWriter {
 
	public:

		GenotypeWriter(const string &, const ushort, const vector<Sample> &, const Chromosomes &, const Filters &);

		void addGenotypes(vector<Genotypes *> *);
		void finalise(const string &, const Chromosomes &, const string &, const OptionsContainer &, const Filters &);

	private:

		const vector<Sample> samples;

		string tmp_filename;
		unordered_map<string, uint> tmp_outfile_chrom_stats;

		ofstream tmp_outfile;
		boost::iostreams::filtering_ostream tmp_outfile_fstream;

    	ProducerConsumerQueue<vector<Genotypes *> * > * genotypes_queue;
		vector<thread> writing_threads; 

		void writeGenotypes(const Chromosomes &, const Filters & filters);

		template <typename ValueType>
		void writeAlleleField(const vector<ValueType> &);

		void writeAlleleSequences(const VariantInfo &, const string &);
		void writeQualityAndFilter(const Genotypes::VariantStats &, const Filters &);
		void writeVariantStats(const Genotypes::VariantStats &, const ushort);
		void writeAlleleCover(vector<ushort> *, const ushort);
		void writeAlleleOrigin(const VariantInfo &);
		void writeSamples(const vector<Genotypes::SampleStats> &, const ushort);
		void writeAlleleKmerStats(const AlleleKmerStats &);

		string generateHeader(const string &, const Chromosomes &, const string &, const string &, const Filters &);

		struct GenotypedVariant {

			uint position;
			string variant_id;
			uint max_ref_length;
			string genotypes;
		};

		static bool genotypedVariantCompare(const GenotypedVariant &, const GenotypedVariant &);
};

#endif