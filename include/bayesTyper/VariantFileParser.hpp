
/*
VariantFileParser.hpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#ifndef __bayesTyper__VariantFileParser_hpp
#define __bayesTyper__VariantFileParser_hpp

#include <string>
#include <unordered_map>
#include <list>
#include <mutex>
#include <condition_variable>
#include <memory>
#include <thread>
#include <set>

#include "boost/iostreams/filtering_stream.hpp"
#include "boost/iostreams/filter/gzip.hpp"

#include "Utils.hpp"
#include "VariantCluster.hpp"
#include "VariantClusterGraph.hpp"
#include "VariantClusterGroup.hpp"
#include "InferenceUnit.hpp"
#include "OptionsContainer.hpp"
#include "Chromosomes.hpp"
#include "ProducerConsumerQueue.hpp"


using namespace std;

class VariantFileParser {
  
	public:

		VariantFileParser(const OptionsContainer &);

		bool constructVariantClusterGroups(InferenceUnit *, const uint, const Chromosomes &);

		struct InterClusterRegion {

			string chrom_name;
			Utils::ChromClass chrom_class;

			uint start;
			uint end;

			InterClusterRegion(const string & chrom_name_in, const Utils::ChromClass & chrom_class_in, const uint start_in, const uint end_in) : chrom_name(chrom_name_in), chrom_class(chrom_class_in), start(start_in), end(end_in) {}
		};

		string getVariantStatsString(const uint);
		
		void sortInterclusterRegions();
		void writeInterclusterRegions(const string &); 

		const vector<InterClusterRegion> & getInterclusterRegions();
		ulong getInterclusterRegionLength();
		ulong getNumberOfInterclusterRegionKmers();
				
		uint getNumberOfVariants();

	private:

		const ushort num_threads;
		const uint max_allele_length;
		const double copy_number_variant_threshold;

		enum class AlleleCount : uchar {Total = 0, Excluded_decoy, Excluded_genome, Excluded_match, Excluded_end, Excluded_length, ALLELE_COUNT_SIZE};

    	vector<uint> allele_type_counter;
    	vector<uint> variant_type_counter;

		uint num_variants;
		uint num_variant_clusters;
	    uint num_variant_cluster_groups;

		vector<InterClusterRegion> intercluster_regions;
		ulong intercluster_regions_length;

	    unordered_set<string> intercluster_chromosomes;

    	string prev_chrom_name;

	    int prev_position;
	    int prev_var_end_position;

		ifstream variants_infile;
		boost::iostreams::filtering_istream variants_infile_fstream;

		uint total_num_variants;
	    
	    bool has_format;
	    vector<string> variant_line;

		void openVariantFile(const string &);
		void updateVariantLine();

		void addSequenceToInterclusterRegions(const string &, const Utils::ChromClass &, const uint, const uint);

		void parseVariants(ProducerConsumerQueue<vector<unordered_map<uint, VariantCluster*> * > * > *, const uint, const Chromosomes &);

		pair<string, bool> getInfoAttributeString(const string &, const string &);
		void rightTrimAllele(string *, string *);

		void addAlternativeAllele(VariantCluster::Variant *, const string &, const string &, const string &);
		VariantCluster::VariantType classifyAllele(const int, const int);		

		uint copyNumberVariantLength(const string &, const string &, const uint);
		void clusterVariants(VariantCluster::Variant &, const uint, const set<uint>, const string &, const Utils::ChromClass &, map<uint, VariantCluster*> *, unordered_map<uint, VariantCluster*> *, list<unordered_set<uint> > *);
		
		void processVariantClusterGroups(ProducerConsumerQueue<vector<unordered_map<uint, VariantCluster*> * > * > *, vector<unordered_map<uint, VariantCluster*> * > **, uint *, unordered_map<uint, VariantCluster*> **, list<unordered_set<uint> > *, map<uint, VariantCluster*> *);
		void mergeVariantClusters(unordered_map<uint, VariantCluster*> *, list<unordered_set<uint> > &);
		
		void processVariantClusterGroupsCallback(vector<VariantClusterGroup*> *, ProducerConsumerQueue<vector<unordered_map<uint, VariantCluster*> * > * > *, mutex *, const Chromosomes &);

		unordered_map<uint, uint> getVariantClusterGroupDependencies(unordered_map<uint, VariantCluster*> *);
};


#endif
