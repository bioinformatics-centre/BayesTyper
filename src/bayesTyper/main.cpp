
/*
main.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#include <time.h>
#include <string>
#include <sstream>
#include <map>
#include <iostream>
#include <assert.h>
#include <fstream>
#include <regex>
#include <iterator>
#include <algorithm>
#include <atomic>
#include <functional>
#include <exception>
#include <bitset>

#include "boost/program_options.hpp"
#include "boost/filesystem.hpp"
#include "boost/serialization/serialization.hpp"
#include "boost/archive/binary_oarchive.hpp"
#include "boost/archive/binary_iarchive.hpp"

#include "boost/iostreams/filtering_streambuf.hpp"
#include "boost/iostreams/filtering_stream.hpp"
#include "boost/iostreams/filter/gzip.hpp"

#include "Utils.hpp"
#include "OptionsContainer.hpp"
#include "KmerHash.hpp"
#include "InferenceUnit.hpp"
#include "VariantClusterGroup.hpp"
#include "InferenceEngine.hpp"
#include "GenotypeWriter.hpp"
#include "KmerCounter.hpp"
#include "Chromosomes.hpp"
#include "Nucleotide.hpp"



using namespace std;
namespace po = boost::program_options;


static const ushort max_num_samples = 30;
static const uint max_parameter_kmers = 1000000;

static const ushort min_filter_samples = 10;

static const string intercluster_regions_file_prefix = "intercluster_regions";
static const string multigroup_kmers_file_prefix = "multigroup_kmers";
static const string parameter_kmers_file_prefix = "parameter_kmers";


int main (int argc, char * const argv[]) {

	assert(Utils::uchar_overflow == (pow(2,8) - 1));
	assert(Utils::ushort_overflow == (pow(2,16) - 1));
	assert(Utils::uint_overflow == (pow(2,32) - 1));
	assert(Utils::ulong_overflow == (pow(2,64) - 1));

	assert(Utils::kmer_size <= static_cast<uint>(Utils::uchar_overflow)); 

	cout << "\n[" << Utils::getLocalTime() << "] " << "You are using BayesTyper (" << BT_VERSION << ")\n" << endl;

	stringstream command_info;

	command_info << "Usage: bayesTyper <command> [options]" << endl;
	command_info << "\nCommands:\n" << endl;
	command_info << "\tcluster\t\tcreate variant clusters" << endl;
	command_info << "\tgenotype\tgenotype variant clusters" << endl;

	if (argc == 1) {

		cout << command_info.str() << endl;

	} else {

		po::options_description help_options("", 160);
		help_options.add_options()
		
		   ("help,h", "produce help message for options")
		;

		if (strcmp(argv[1], "cluster") == 0) {

			OptionsContainer options_container("cluster", BT_VERSION, Utils::getLocalTime());

			po::options_description required_options("== Required ==", 160);
			required_options.add_options()

				("variant-file,v", po::value<string>()->required()->notifier(bind(&OptionsContainer::parseValue<string>, &options_container, "variant-file", placeholders::_1)), "variant file (vcf format).")
			 	("samples-file,s", po::value<string>()->required()->notifier(bind(&OptionsContainer::parseValue<string>, &options_container, "samples-file", placeholders::_1)), "samples file (see github documentation for format specifications).")
				("genome-file,g", po::value<string>()->required()->notifier(bind(&OptionsContainer::parseValue<string>, &options_container, "genome-file", placeholders::_1)), "reference genome file (fasta format).")
			;

			po::options_description general_options("== General ==", 160);
			general_options.add_options()

				("decoy-file,d", po::value<string>()->default_value("")->notifier(bind(&OptionsContainer::parseValue<string>, &options_container, "decoy-file", placeholders::_1)), "decoy sequences file (fasta format).")
				("output-prefix,o", po::value<string>()->default_value("bayestyper")->notifier(bind(&OptionsContainer::parseValue<string>, &options_container, "output-prefix", placeholders::_1)), "output prefix.")
		    	("random-seed,r", po::value<uint>()->default_value(time(nullptr), "unix time")->notifier(bind(&OptionsContainer::parseValue<uint>, &options_container, "random-seed", placeholders::_1)), "seed for pseudo-random number generator.")
				("threads,p", po::value<ushort>()->default_value(1)->notifier(bind(&OptionsContainer::parseValue<ushort>, &options_container, "threads", placeholders::_1)), "number of threads used (+= 2 I/O threads).")
			;

			po::options_description parameters_options("== Parameters ==", 160);
			parameters_options.add_options()

				("min-number-of-unit-variants", po::value<uint>()->default_value(5000000)->notifier(bind(&OptionsContainer::parseValue<uint>, &options_container, "min-number-of-unit-variants", placeholders::_1)), "minimum number of variants per inference unit.")
				("max-allele-length", po::value<uint>()->default_value(500000)->notifier(bind(&OptionsContainer::parseValue<uint>, &options_container, "max-allele-length", placeholders::_1)), "exclude alleles (reference and alternative) longer than <length>.")
				("copy-number-variant-threshold", po::value<float>()->default_value(0.5, "0.5")->notifier(bind(&OptionsContainer::parseValue<float>, &options_container, "copy-number-variant-threshold", placeholders::_1)), "minimum fraction of identical kmers required between an allele and the downstream reference sequence in order for it to be classified as a copy number.")
				("max-number-of-sample-haplotype-candidates", po::value<ushort>()->default_value(16)->notifier(bind(&OptionsContainer::parseValue<ushort>, &options_container, "max-number-of-sample-haplotype-candidates", placeholders::_1)), "maximum number of haplotype candidates per sample.")
			;

			po::options_description desc("## BayesTyper cluster options ##");
			desc.add(help_options).add(required_options).add(general_options).add(parameters_options);

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);

			if (vm.count("help") or (argc == 2)) {
				
				cout << desc << endl;
				return 1;
			}

			po::notify(vm);

		    assert(options_container.getValue<uint>("min-number-of-unit-variants") > 0);

		    assert(options_container.getValue<float>("copy-number-variant-threshold") >= 0);
		    assert(options_container.getValue<float>("copy-number-variant-threshold") <= 1);
		    assert(options_container.getValue<ushort>("max-number-of-sample-haplotype-candidates") > 0);
		    assert((options_container.getValue<ushort>("max-number-of-sample-haplotype-candidates") * max_num_samples) < Utils::ushort_overflow);


			cout << "[" << Utils::getLocalTime() << "] Seeding pseudo-random number generator with " << options_container.getValue<uint>("random-seed") << " ..." << endl;
			cout << "[" << Utils::getLocalTime() << "] Setting the kmer size to " << Utils::kmer_size << " ..." << endl;

			vector<Sample> samples;

			ifstream samples_infile(options_container.getValue<string>("samples-file").c_str());

	        if (!samples_infile.is_open()) {

	            cerr << "\nERROR: Unable to open file " << options_container.getValue<string>("samples-file").c_str() << "\n" << endl;
	            exit(1);
	        }

			for (string sample_line; getline(samples_infile, sample_line);) {

				samples.emplace_back(sample_line);
			}

		    samples_infile.close();

		    if (samples.empty()) {

				cerr << "\nERROR: Samples file empty\n" << endl;
				exit(1);    	
		    }

		    if (samples.size() > max_num_samples) {

				cerr << "\nERROR: The maximum number of samples supported by BayesTyper is currently " << max_num_samples << "\n" << endl;
				exit(1);    	
		    }

		    const ushort num_threads = options_container.getValue<ushort>("threads");
		    const string output_prefix = options_container.getValue<string>("output-prefix");

			cout << "\n[" << Utils::getLocalTime() << "] Parsed information for " << samples.size() << " sample(s)" << endl;
		    cout << "\n[" << Utils::getLocalTime() << "] Parsing reference genome ..." << endl;

		    Chromosomes chromosomes(options_container.getValue<string>("genome-file"), false);

		    cout << "[" << Utils::getLocalTime() << "] Parsed " << chromosomes.getTotalCount() << " reference genome chromosomes(s) (" << chromosomes.getTotalLength() << " nucleotides)" << endl;

		    cout << "\n[" << Utils::getLocalTime() << "] Parsing decoy sequence(s) ..." << endl;

		    chromosomes.addFasta(options_container.getValue<string>("decoy-file"), true);

		    cout << "[" << Utils::getLocalTime() << "] Parsed " << chromosomes.getDecoyCount() << " decoy sequence(s) (" << chromosomes.getDecoyLength() << " nucleotides)" << endl;

		    chromosomes.convertToUpper();
		    assert(chromosomes.getTotalLength() >= chromosomes.getDecoyLength());


			KmerCounter kmer_counter(samples, options_container);
			VariantFileParser variant_file_parser(options_container);
		    
		    const uint num_variants = variant_file_parser.getNumberOfVariants();
		    const uint num_units = max(uint(1), static_cast<uint>(floor(num_variants/static_cast<float>(options_container.getValue<uint>("min-number-of-unit-variants")))));

		    cout << "\n[" << Utils::getLocalTime() << "] Setting the number of inference units to " << num_units << " across " << num_variants << " variants ..." << endl;

			const ulong expected_num_path_kmers = ceil((chromosomes.getTotalLength() - chromosomes.getDecoyLength()) * (1 + (0.05 * 2 * samples.size())));
			ulong num_path_kmers = 0;
			
			ThreadedKmerBloom<Utils::kmer_size> path_kmer_bloom(expected_num_path_kmers, 0.0001);
			BooleanKmerHash multigroup_kmer_hash(ceil(expected_num_path_kmers * 0.01), num_threads);

			cout << "\n[" << Utils::getLocalTime() << "] " << Utils::getMaxMemoryUsage() << endl;


			bool variant_file_parsed = false;

		    for (uint unit_idx = 1; unit_idx < (num_units + 1); unit_idx++) {

		      	cout << "\n" << endl;

		    	assert(!variant_file_parsed);

				InferenceUnit inference_unit(unit_idx);
				inference_unit.cluster_options_header = options_container.getHeader();

		    	variant_file_parsed = variant_file_parser.constructVariantClusterGroups(&inference_unit, ceil(num_variants/static_cast<float>(num_units)), chromosomes);

		    	assert(inference_unit.num_variant_clusters <= inference_unit.num_variants);
		    	assert(inference_unit.variant_cluster_groups.size() <= inference_unit.num_variant_clusters);

				sort(inference_unit.variant_cluster_groups.begin(), inference_unit.variant_cluster_groups.end(), VariantClusterGroupCompare);

		      	cout << endl;

		    	kmer_counter.findVariantClusterPaths(&inference_unit, options_container.getValue<ushort>("max-number-of-sample-haplotype-candidates"));
		    	kmer_counter.countPathMultigroupKmers(&multigroup_kmer_hash, &path_kmer_bloom, &inference_unit);

				num_path_kmers += inference_unit.num_path_kmers;

		   		boost::filesystem::path inference_unit_dir(output_prefix + "_unit_" + to_string(inference_unit.index));
		  		
		  		if (!boost::filesystem::create_directory(inference_unit_dir)) {

		  			cerr << "\nERROR: Unit directory " << inference_unit_dir.string() << "/ already exist\n" << endl;
		  			exit(1);
		  		}

		  		string clusters_filename = inference_unit_dir.string() + "/variant_clusters.bin";

				{
					ofstream clusters_outfile(clusters_filename, ios::binary);

				    if (!clusters_outfile.is_open()) {

				        cerr << "\nERROR: Unable to write file " << clusters_filename << "\n" << endl;
				        exit(1);
				    }

			  		boost::iostreams::filtering_ostream clusters_outfile_fstream;
					
					clusters_outfile_fstream.push(boost::iostreams::gzip_compressor());
					clusters_outfile_fstream.push(boost::ref(clusters_outfile));

   			 		assert(clusters_outfile_fstream.is_complete());

			        boost::archive::binary_oarchive clusters_outfile_archive(clusters_outfile_fstream);
			        clusters_outfile_archive << inference_unit;
			    }

			    cout << "\n[" << Utils::getLocalTime() << "] Wrote unit " << inference_unit.index << " variant clusters to " << clusters_filename << endl;

				for (auto & variant_cluster_group: inference_unit.variant_cluster_groups) {

					delete variant_cluster_group;
				}

				cout << "\n[" << Utils::getLocalTime() << "] " << Utils::getMaxMemoryUsage() << endl;
		    }

		    assert(variant_file_parsed);


		    cout << "\n\n" << variant_file_parser.getVariantStatsString(num_units);
    
		    if (expected_num_path_kmers < num_path_kmers) {

		    	cout << "\nWARNING: Multigroup kmer estimate might be inflated due to the number of kmers being higher than expected.\n" << endl;
		    }
		    
	   		boost::filesystem::path cluster_data_dir(output_prefix + "_cluster_data");
	  		
	  		if (!boost::filesystem::create_directory(cluster_data_dir)) {

	  			cerr << "\nERROR: Cluster data directory " << cluster_data_dir.string() << "/ already exist\n" << endl;
	  			exit(1);
	  		}

		    cout << "\n[" << Utils::getLocalTime() << "] " << Utils::getMaxMemoryUsage() << endl;
    		

    		cout << "\n\n[" << Utils::getLocalTime() << "] Writing inter-cluster regions ..." << endl;

    		const string intercluster_regions_dir_prefix = cluster_data_dir.string() + "/" + intercluster_regions_file_prefix;

		    variant_file_parser.sortInterclusterRegions();
		    variant_file_parser.writeInterclusterRegions(intercluster_regions_dir_prefix);

    		cout << "[" << Utils::getLocalTime() << "] Wrote " << variant_file_parser.getInterclusterRegions().size() << " regions to " << intercluster_regions_dir_prefix << ".txt.gz\n" << endl;


	      	const uint max_intercluster_kmers = 3 * max_parameter_kmers;
	      	const float parameter_kmer_fraction = min(float(1), static_cast<float>(max_intercluster_kmers) / variant_file_parser.getNumberOfInterclusterRegionKmers());
	
    		BooleanKmerHash parameter_kmer_hash(max_intercluster_kmers + chromosomes.getDecoyLength(), num_threads);

			kmer_counter.countInterclusterParameterKmers(&parameter_kmer_hash, variant_file_parser.getInterclusterRegions(), chromosomes, path_kmer_bloom, parameter_kmer_fraction);

		    parameter_kmer_hash.shuffle(options_container.getValue<uint>("random-seed"));

		    const string parameter_kmers_dir_prefix = cluster_data_dir.string() + "/" + parameter_kmers_file_prefix;

		    auto num_parameter_kmers = parameter_kmer_hash.writeKmersToFasta(parameter_kmers_dir_prefix, true, max_parameter_kmers);
		    assert(num_parameter_kmers <= max_parameter_kmers);

		    cout << "[" << Utils::getLocalTime() << "] Wrote " << num_parameter_kmers << " kmers to " << parameter_kmers_dir_prefix << ".fa.gz" << endl;


			cout << "\n[" << Utils::getLocalTime() << "] Creating multigroup kmers bloom filter ..." << endl;
	    
			const ulong num_multigroup_kmers = multigroup_kmer_hash.size();
		    KmerBloom<Utils::kmer_size> multigroup_kmer_bloom(num_multigroup_kmers, 0.0001);

		    assert(num_multigroup_kmers == multigroup_kmer_hash.addKmersToBloomFilter(&multigroup_kmer_bloom, false));

		    const string multigroup_kmers_dir_prefix = cluster_data_dir.string() + "/" + multigroup_kmers_file_prefix;
		    multigroup_kmer_bloom.save(multigroup_kmers_dir_prefix);

		    cout << "[" << Utils::getLocalTime() << "] Wrote " << num_multigroup_kmers << " kmers to " << multigroup_kmers_dir_prefix << ".bloom[Meta|Data]" << endl;

			cout << "\n[" << Utils::getLocalTime() << "] " << Utils::getMaxMemoryUsage() << endl;
			

			cout << "\n\n[" << Utils::getLocalTime() << "] BayesTyper cluster completed succesfully!\n" << endl;
		
		} else if (strcmp(argv[1], "genotype") == 0) {

			const string cluster_data_dir_option_help = "cluster data directory containing " + intercluster_regions_file_prefix + ".txt.gz, " + multigroup_kmers_file_prefix + ".bloom[Meta|Data] & " + parameter_kmers_file_prefix + ".fa.gz (BayesTyper cluster output).";

			const string min_homozygote_genotypes_option_help = "filter variants with less than <value> homozygote genotypes (calculated before other filters). Minimum " + to_string(min_filter_samples) + " samples required for this filter.";


			OptionsContainer options_container("genotype", BT_VERSION, Utils::getLocalTime());

			po::options_description required_options("== Required ==", 160);
			required_options.add_options()

				("variant-clusters-file,v", po::value<string>()->required()->notifier(bind(&OptionsContainer::parseValue<string>, &options_container, "variant-clusters-file", placeholders::_1)), "variant_clusters.bin file (BayesTyper cluster output).")
				("cluster-data-dir,c", po::value<string>()->required()->notifier(bind(&OptionsContainer::parseValue<string>, &options_container, "cluster-data-dir", placeholders::_1)), cluster_data_dir_option_help.c_str())
			 	("samples-file,s", po::value<string>()->required()->notifier(bind(&OptionsContainer::parseValue<string>, &options_container, "samples-file", placeholders::_1)), "samples file (see github documentation for format specifications).")
				("genome-file,g", po::value<string>()->required()->notifier(bind(&OptionsContainer::parseValue<string>, &options_container, "genome-file", placeholders::_1)), "reference genome file (fasta format).")
			;

			po::options_description general_options("== General ==", 160);
			general_options.add_options()

				("decoy-file,d", po::value<string>()->default_value("")->notifier(bind(&OptionsContainer::parseValue<string>, &options_container, "decoy-file", placeholders::_1)), "decoy sequences file (fasta format).")
				("output-prefix,o", po::value<string>()->default_value("bayestyper")->notifier(bind(&OptionsContainer::parseValue<string>, &options_container, "output-prefix", placeholders::_1)), "output prefix.")
				("gzip-output,z", po::value<bool>()->default_value(false)->implicit_value(true)->notifier(bind(&OptionsContainer::parseValue<bool>, &options_container, "gzip-output", placeholders::_1)), "compress <output-prefix>.vcf using gzip.")
		    	("random-seed,r", po::value<uint>()->default_value(time(nullptr), "unix time")->notifier(bind(&OptionsContainer::parseValue<uint>, &options_container, "random-seed", placeholders::_1)), "seed for pseudo-random number generator.")
				("threads,p", po::value<ushort>()->default_value(1)->notifier(bind(&OptionsContainer::parseValue<ushort>, &options_container, "threads", placeholders::_1)), "number of threads used (+= 2 I/O threads).")
				("chromosome-ploidy-file,y", po::value<string>()->default_value("")->notifier(bind(&OptionsContainer::parseValue<string>, &options_container, "chromosome-ploidy-file", placeholders::_1)), "chromosome gender ploidy file (see github documentation for format specifications). Human ploidy levels will be assumed if no file is given.")
			;

			po::options_description genotyping_parameters_options("== Genotyping parameters ==", 160);
			genotyping_parameters_options.add_options()

				("gibbs-burn-in", po::value<ushort>()->default_value(100)->notifier(bind(&OptionsContainer::parseValue<ushort>, &options_container, "gibbs-burn-in", placeholders::_1)), "number of burn-in iterations.")
				("gibbs-samples", po::value<ushort>()->default_value(250)->notifier(bind(&OptionsContainer::parseValue<ushort>, &options_container, "gibbs-samples", placeholders::_1)), "number of Gibbs iterations.")
				("number-of-gibbs-chains", po::value<ushort>()->default_value(20)->notifier(bind(&OptionsContainer::parseValue<ushort>, &options_container, "number-of-gibbs-chains", placeholders::_1)), "number of parallel Gibbs sampling chains.")
				("kmer-subsampling-rate", po::value<float>()->default_value(0.1, "0.1")->notifier(bind(&OptionsContainer::parseValue<float>, &options_container, "kmer-subsampling-rate", placeholders::_1)), "subsampling rate for subsetting kmers used for genotype inference (a new subset is sampled for each Gibbs sampling chain).")
				("max-haplotype-variant-kmers", po::value<uint>()->default_value(500)->notifier(bind(&OptionsContainer::parseValue<uint>, &options_container, "max-haplotype-variant-kmers", placeholders::_1)), "maximum number of kmers used for genotype inference after subsampling across a haplotype candidate for each variant (a new subset is sampled for each Gibbs sampling chain).")
				("noise-rate-prior", po::value<string>()->default_value("1,1")->notifier(bind(&OptionsContainer::parseValuePair<float>, &options_container, "noise-rate-prior", placeholders::_1)), "parameters for Poisson noise rate gamma prior (<shape>,<scale>). All samples will use the same parameters.")
			;

			po::options_description filter_parameters_options("== Filter parameters ==", 160);
			filter_parameters_options.add_options()

				("min-homozygote-genotypes", po::value<ushort>()->default_value(1)->notifier(bind(&OptionsContainer::parseValue<ushort>, &options_container, "min-homozygote-genotypes", placeholders::_1)), min_homozygote_genotypes_option_help.c_str())
				("min-genotype-posterior", po::value<float>()->default_value(0.99, "0.99")->notifier(bind(&OptionsContainer::parseValue<float>, &options_container, "min-genotype-posterior", placeholders::_1)), "filter genotypes with a posterior probability (GPP) below <value>.")
				("min-number-of-kmers", po::value<float>()->default_value(1, "1")->notifier(bind(&OptionsContainer::parseValue<float>, &options_container, "min-number-of-kmers", placeholders::_1)), "filter sampled alleles with less than <value> kmers (NAK).")
				("disable-observed-kmers", po::value<bool>()->default_value(false)->implicit_value(true)->notifier(bind(&OptionsContainer::parseValue<bool>, &options_container, "disable-observed-kmers", placeholders::_1)), "disable filtering of sampled alleles with a low fraction of observed kmers (FAK).")
			;

			po::options_description desc("## BayesTyper genotype options ##");
			desc.add(help_options).add(required_options).add(general_options).add(genotyping_parameters_options).add(filter_parameters_options);

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);

			if (vm.count("help") or (argc == 2)) {
				
				cout << desc << endl;
				return 1;
			}

			po::notify(vm);

		    assert(options_container.getValue<ushort>("gibbs-burn-in") > 0);
		    assert(options_container.getValue<ushort>("gibbs-samples") > 0);
		    assert(options_container.getValue<ushort>("number-of-gibbs-chains") > 0);
		    assert(options_container.getValue<float>("kmer-subsampling-rate") > 0);
		    assert(options_container.getValue<float>("kmer-subsampling-rate") <= 1);

		    assert((options_container.getValue<pair<float,float> >("noise-rate-prior").first) > 0);
		    assert((options_container.getValue<pair<float,float> >("noise-rate-prior").second) > 0);


			cout << "[" << Utils::getLocalTime() << "] Seeding pseudo-random number generator with " << options_container.getValue<uint>("random-seed") << " ..." << endl;
			cout << "[" << Utils::getLocalTime() << "] Setting the kmer size to " << Utils::kmer_size << " ..." << endl;

			vector<Sample> samples;

			ifstream samples_infile(options_container.getValue<string>("samples-file").c_str());

	        if (!samples_infile.is_open()) {

	            cerr << "\nERROR: Unable to open file " << options_container.getValue<string>("samples-file").c_str() << "\n" << endl;
	            exit(1);
	        }

			for (string sample_line; getline(samples_infile, sample_line);) {

				samples.emplace_back(sample_line);
			}

		    samples_infile.close();

		    if (samples.empty()) {

				cerr << "\nERROR: Samples file empty\n" << endl;
				exit(1);    	
		    }

		    if (samples.size() > max_num_samples) {

				cerr << "\nERROR: The maximum number of samples supported by BayesTyper is currently " << max_num_samples << "\n" << endl;
				exit(1);    	
		    }

		    const ushort num_threads = options_container.getValue<ushort>("threads");
		    const string output_prefix = options_container.getValue<string>("output-prefix");

			cout << "\n[" << Utils::getLocalTime() << "] Parsed information for " << samples.size() << " sample(s)" << endl;

			KmerCounter kmer_counter(samples, options_container);

		    cout << "\n[" << Utils::getLocalTime() << "] Parsing reference genome ..." << endl;

		    Chromosomes chromosomes(options_container.getValue<string>("genome-file"), false);

		    cout << "[" << Utils::getLocalTime() << "] Parsed " << chromosomes.getTotalCount() << " reference genome chromosomes(s) (" << chromosomes.getTotalLength() << " nucleotides)" << endl;

		    cout << "\n[" << Utils::getLocalTime() << "] Parsing decoy sequence(s) ..." << endl;

		    chromosomes.addFasta(options_container.getValue<string>("decoy-file"), true);

		    cout << "[" << Utils::getLocalTime() << "] Parsed " << chromosomes.getDecoyCount() << " decoy sequence(s) (" << chromosomes.getDecoyLength() << " nucleotides)" << endl;

		    chromosomes.convertToUpper();
		    assert(chromosomes.getTotalLength() >= chromosomes.getDecoyLength());

	    	cout << "\n[" << Utils::getLocalTime() << "] " << Utils::getMaxMemoryUsage() << endl;
		    

		    cout << "\n\n[" << Utils::getLocalTime() << "] Parsing variant clusters ..." << endl;

			InferenceUnit inference_unit;

			{
				std::ifstream clusters_infile(options_container.getValue<string>("variant-clusters-file"), std::ios::binary);

		        if (!clusters_infile.is_open()) {

		            cerr << "\nERROR: Unable to open file " << options_container.getValue<string>("variant-clusters-file") << "\n" << endl;
		            exit(1);
		        }

			  	boost::iostreams::filtering_istream clusters_infile_fstream;
				
				clusters_infile_fstream.push(boost::iostreams::gzip_decompressor());
	    		clusters_infile_fstream.push(boost::ref(clusters_infile));

    			assert(clusters_infile_fstream.is_complete());    

		        boost::archive::binary_iarchive clusters_infile_archive(clusters_infile_fstream);
		        clusters_infile_archive >> inference_unit;   
			}

		    cout << "[" << Utils::getLocalTime() << "] Parsed " << inference_unit.num_variant_clusters << " variant clusters (" << inference_unit.num_variants << " variants)" << endl;

    		const string intercluster_regions_dir_prefix = options_container.getValue<string>("cluster-data-dir") + "/" + intercluster_regions_file_prefix;
		    const string parameter_kmers_dir_prefix = options_container.getValue<string>("cluster-data-dir") + "/" + parameter_kmers_file_prefix;
			const string multigroup_kmers_dir_prefix = options_container.getValue<string>("cluster-data-dir") + "/" + multigroup_kmers_file_prefix;

		    ThreadedKmerBloom<Utils::kmer_size> * path_kmer_bloom = new ThreadedKmerBloom<Utils::kmer_size>(inference_unit.num_path_kmers + max_parameter_kmers, 0.0001);

		    KmerCountsHash * kmer_hash;

			if (samples.size() < 4) {

		        kmer_hash = new ObservedKmerCountsHash<3>(inference_unit.num_path_kmers + max_parameter_kmers, num_threads); 

			} else if (samples.size() < 11) {

		        kmer_hash = new ObservedKmerCountsHash<10>(inference_unit.num_path_kmers + max_parameter_kmers, num_threads); 

		    } else if ((samples.size() > 10) and (samples.size() < 21)) {

		        kmer_hash = new ObservedKmerCountsHash<20>(inference_unit.num_path_kmers + max_parameter_kmers, num_threads); 

		    } else {

		        assert(samples.size() <= 30);
		        kmer_hash = new ObservedKmerCountsHash<30>(inference_unit.num_path_kmers + max_parameter_kmers, num_threads);
		    }

		    cout << "\n[" << Utils::getLocalTime() << "] Parsing parameter kmers ..." << endl;

		    uint num_parameter_kmers = 0;

		    {
			    ifstream kmers_infile(parameter_kmers_dir_prefix + ".fa.gz", std::ios::binary);

		        if (!kmers_infile.is_open()) {

		            cerr << "\nERROR: Unable to open file " << parameter_kmers_dir_prefix + ".fa.gz" << "\n" << endl;
		            exit(1);
		        }

				boost::iostreams::filtering_istream kmers_infile_fstream;
				
				kmers_infile_fstream.push(boost::iostreams::gzip_decompressor());
		    	kmers_infile_fstream.push(boost::ref(kmers_infile));

    			assert(kmers_infile_fstream.is_complete());   

			    string parameter_kmer_str;
			    getline(kmers_infile_fstream, parameter_kmer_str);

			    assert(parameter_kmer_str == (">k" + to_string(Utils::kmer_size)));

			    while (getline(kmers_infile_fstream, parameter_kmer_str)) {

			    	num_parameter_kmers++;
			    	assert(parameter_kmer_str.size() == Utils::kmer_size);

			    	path_kmer_bloom->addKmer(parameter_kmer_str);

			    	auto parameter_kmer = Nucleotide::ntToBit<Utils::kmer_size>(parameter_kmer_str);
			    	assert(parameter_kmer.second);

					auto kmer_counts = kmer_hash->addKmer(parameter_kmer.first, false);

					assert(kmer_counts.first);
					assert(kmer_counts.second);

					kmer_counts.first->isParameter(true);
				}
			}

		    assert(num_parameter_kmers <= max_parameter_kmers);
			kmer_hash->sortKmers();

		    cout << "[" << Utils::getLocalTime() << "] Parsed " << num_parameter_kmers << " kmers" << endl;
	      	cout << "\n[" << Utils::getLocalTime() << "] " << Utils::getMaxMemoryUsage() << endl;
	      

	      	cout << "\n" << endl;

	      	auto chrom_ploidy = ChromosomePloidy(options_container.getValue<string>("chromosome-ploidy-file"), chromosomes, samples);

		    kmer_counter.countPathKmers(path_kmer_bloom, &inference_unit);
	    	kmer_counter.countInterclusterKmers(kmer_hash, path_kmer_bloom, intercluster_regions_dir_prefix, chromosomes, chrom_ploidy);

	    	cout << endl;
	      	kmer_counter.parseSampleKmers(kmer_hash, path_kmer_bloom);

	      	delete path_kmer_bloom;
		
			cout << endl;	
			kmer_counter.classifyPathKmers(kmer_hash, &inference_unit, multigroup_kmers_dir_prefix);
			
			auto intercluster_kmer_stats = kmer_hash->calculateKmerStats(samples);

			cout << "\n[" << Utils::getLocalTime() << "] " << Utils::getMaxMemoryUsage() << endl;
	      	

	      	cout << "\n" << endl;

			CountDistribution count_distribution(samples, options_container);
			count_distribution.setGenomicCountDistributions(intercluster_kmer_stats, output_prefix + "_genomic_parameters");

			cout << "\n[" << Utils::getLocalTime() << "] " << Utils::getMaxMemoryUsage() << endl;
	      	

	      	cout << "\n" << endl;

			InferenceEngine inference_engine(samples, chrom_ploidy, options_container);
			inference_engine.estimateNoiseParameters(&count_distribution, &inference_unit, kmer_hash, output_prefix + "_noise_parameters");	

			cout << "\n[" << Utils::getLocalTime() << "] " << Utils::getMaxMemoryUsage() << endl;
	      	

	      	cout << "\n" << endl;

			Filters filters(options_container, count_distribution.getGenomicCountDistributions(), min_filter_samples);
			GenotypeWriter genotype_writer(output_prefix, num_threads, samples, chromosomes, filters);

			inference_engine.genotypeVariantClusterGroups(&inference_unit, kmer_hash, count_distribution, filters, &genotype_writer);

			delete kmer_hash;

			cout << "\n[" << Utils::getLocalTime() << "] " << Utils::getMaxMemoryUsage() << endl;


	      	cout << "\n" << endl;

			genotype_writer.finalise(output_prefix, chromosomes, inference_unit.cluster_options_header, options_container, filters);		

			cout << "\n[" << Utils::getLocalTime() << "] " << Utils::getMaxMemoryUsage() << endl;


			cout << "\n\n[" << Utils::getLocalTime() << "] BayesTyper genotype completed succesfully!\n" << endl;

		} else {

			cout << command_info.str() << endl;
		}
	}

	return 0;
}
