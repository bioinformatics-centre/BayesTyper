
/*
split.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#include <sstream>
#include <math.h>
#include <thread>

#include "boost/filesystem.hpp"

#include "JoiningString.hpp"
#include "Auxiliaries.hpp"

#include "split.hpp"

namespace Split {

	struct Batch {

		uint start_variant_idx;
		string directory;
	};

	void splitCallback(const string & vcf_batch_name, const string & vcf_file, const vector<Batch> & batches) {

		GenotypedVcfFileReader vcf_reader(vcf_file, true);
		auto output_meta_data = vcf_reader.metaData();

		auto input_filename = Utils::splitString(vcf_file, '/').back();
		assert(input_filename.size() > 4);
		assert(input_filename.substr(input_filename.size() - 4, 4) == ".vcf");

		const string output_filename = input_filename.substr(0, input_filename.size() - 4) + "_" + vcf_batch_name + ".vcf";

		uint num_vars = 0;

		Variant * cur_var;
		VcfFileWriter * vcf_writer = nullptr;

		auto batch_it = batches.begin();
		
		assert(batch_it != batches.end());
		assert(batch_it->start_variant_idx == 0);

		string cur_variant_cluster_group_id = "";
		string prev_variant_cluster_group_id = "";

		while (vcf_reader.getNextVariant(&cur_var)) {
				
			auto vcgi_value = cur_var->info().getValue<string>("VCGI");
			assert(vcgi_value.second);

			cur_variant_cluster_group_id = vcgi_value.first;

			if (batch_it != batches.end()) {

				if (num_vars == batch_it->start_variant_idx) {

					assert(cur_variant_cluster_group_id != prev_variant_cluster_group_id);

					if (batch_it != batches.begin()) {

						delete vcf_writer;
					}

					vcf_writer = new VcfFileWriter(batch_it->directory + "/" + output_filename, output_meta_data, true);
					batch_it++;
				} 
			}

			assert(vcf_writer);
			vcf_writer->write(cur_var);

			num_vars++;
			prev_variant_cluster_group_id = cur_variant_cluster_group_id;

			delete cur_var;
		}

		assert(batch_it == batches.end());

		delete vcf_writer;
	}

	void split(const vector<string> & vcf_filenames, const uint min_batch_size) {

		cout << "[" << Utils::getLocalTime() << "] Running BayesTyperTools (" << BT_VERSION << ") split ...\n" << endl;

		assert(!(vcf_filenames.empty()));

		auto tmpl_vcf_split = Utils::splitString(vcf_filenames.front(), ':');	
		assert(tmpl_vcf_split.size() == 2);

		VcfFileReader tmpl_vcf_reader(tmpl_vcf_split.back(), true);
		Auxiliaries::removeNonRelevantInfoDescriptors(&(tmpl_vcf_reader.metaData()), {"VCGI"});

		cout << "[" << Utils::getLocalTime() << "] Estimating variant batches ..." << endl;

		Variant * cur_var;

		uint num_vars = 0;
		uint cur_batch_size = min_batch_size;

		string cur_variant_cluster_group_id = "";
		string prev_variant_cluster_group_id = "";

		vector<Batch> batches;

		while (tmpl_vcf_reader.getNextVariant(&cur_var)) {
				
			auto vcgi_value = cur_var->info().getValue<string>("VCGI");
			assert(vcgi_value.second);

			cur_variant_cluster_group_id = vcgi_value.first;

			if (cur_variant_cluster_group_id != prev_variant_cluster_group_id) {

				if (cur_batch_size >= min_batch_size) {
				
					batches.emplace_back(Batch());
					batches.back().start_variant_idx = num_vars;

					cur_batch_size = 0;
				} 
			}

			num_vars++;
			cur_batch_size++;

			prev_variant_cluster_group_id = cur_variant_cluster_group_id;

			delete cur_var;
		}

		for (uint batch_idx = 0; batch_idx < batches.size(); batch_idx++) {

			boost::filesystem::path cur_batch_dir("./var_batch_" + to_string(batch_idx + 1));
			assert(boost::filesystem::create_directory(cur_batch_dir));

			batches.at(batch_idx).directory = cur_batch_dir.string();
		}

		cout << "[" << Utils::getLocalTime() << "] Splitting " << num_vars << " variants in " << vcf_filenames.size() << " file(s) into " << batches.size() << " variant batches ..." << endl;

		vector<thread> splitting_threads;

		for (auto &vcf: vcf_filenames) {

			auto vcf_split = Utils::splitString(vcf, ':');	
			assert(vcf_split.size() == 2);

	   	    splitting_threads.emplace_back(thread(&splitCallback, vcf_split.front(), vcf_split.back(), batches));
	    }

	    for (auto & thread: splitting_threads) {
	        	
	       	thread.join();
		}

		cout << "[" << Utils::getLocalTime() << "] Completed splitting\n" << endl;
	}
}
