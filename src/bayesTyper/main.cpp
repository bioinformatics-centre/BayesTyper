
/*
main.cpp - This file is part of BayesTyper (v1.1)


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
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string.hpp"

#include "Utils.hpp"
#include "OptionsContainer.hpp"
#include "KmerHash.hpp"
#include "VariantClusterGroup.hpp"
#include "InferenceEngine.hpp"
#include "GenotypeWriter.hpp"
#include "KmerFactory.hpp"

using namespace std;
namespace po = boost::program_options;

static const ushort default_kmer_size = 55;

int main (int argc, char * const argv[]) {

	OptionsContainer options_container(BT_VERSION, Utils::getLocalTime());

	cout << "\n[" << Utils::getLocalTime() << "] " << "You are using BayesTyper (" << BT_VERSION << ")\n" << endl;
	po::options_description help("", 160);

	help.add_options()
	   ("help,h", "produce help message for options")
	;

	po::options_description required("== Required ==", 160);
	required.add_options()

		("vcf-file,v", po::value<string>()->required()->notifier(bind(&OptionsContainer::parseValue<string>, &options_container, "vcf-file", placeholders::_1)), "variant file.")
	 	("samples-file,s", po::value<string>()->required()->notifier(bind(&OptionsContainer::parseValue<string>, &options_container, "samples-file", placeholders::_1)), "samples file.")
		("genome-file,g", po::value<string>()->required()->notifier(bind(&OptionsContainer::parseValue<string>, &options_container, "genome-file", placeholders::_1)), "reference genome.")
	;

	po::options_description general("== General ==", 160);
	general.add_options()

		("output-prefix,o", po::value<string>()->default_value("bayestyper")->notifier(bind(&OptionsContainer::parseValue<string>, &options_container, "output-prefix", placeholders::_1)), "prefix for output files.")
		("decoy-file,d", po::value<string>()->default_value("")->notifier(bind(&OptionsContainer::parseValue<string>, &options_container, "decoy-file", placeholders::_1)), "decoy sequences.")
    	("random-seed,r", po::value<uint>()->default_value(time(nullptr), "unix time")->notifier(bind(&OptionsContainer::parseValue<uint>, &options_container, "random-seed", placeholders::_1)), "seed for pseudo-random number generator.")
		("threads,p", po::value<ushort>()->default_value(1)->notifier(bind(&OptionsContainer::parseValue<ushort>, &options_container, "threads", placeholders::_1)), "number of threads used (+= 2 I/O threads).")
		("chromosome-regions,c", po::value<string>()->default_value("")->notifier(bind(&OptionsContainer::parseRegions, &options_container, "chromosome-regions", placeholders::_1)), "genotype variants only in the supplied region list (<chr1>:<start1>-<end1>,<chr2>,...); for whole chromosomes the start and end positions can be omitted. If no list is given all variant will be genotyped.")
	;

	po::options_description graph("== Graph ==", 160);
	graph.add_options()

		("kmer-size", po::value<ushort>()->default_value(default_kmer_size)->notifier(bind(&OptionsContainer::parseValue<ushort>, &options_container, "kmer-size", placeholders::_1)), "kmer size (31, 39, 47, 55 or 63).")
		("max-allele-length", po::value<uint>()->default_value(3000000)->notifier(bind(&OptionsContainer::parseValue<uint>, &options_container, "max-allele-length", placeholders::_1)), "exclude alleles (reference and alternative) longer than <length>.")
		("copy-number-variant-threshold", po::value<double>()->default_value(0.5, "0.5")->notifier(bind(&OptionsContainer::parseValue<double>, &options_container, "copy-number-variant-threshold", placeholders::_1)), "minimum fraction of identical kmers required between an allele and the downstream reference sequence in order for it to be classified as a copy number")
		("max-number-of-sample-haplotype-candidates", po::value<ushort>()->default_value(24)->notifier(bind(&OptionsContainer::parseValue<ushort>, &options_container, "max-number-of-sample-haplotype-candidates", placeholders::_1)), "maximum number of haplotype candidates per sample (the total number of candidates will never exceed 256).")
	;

	po::options_description inference("== Inference ==", 160);
	inference.add_options()

		("gibbs-burn-in", po::value<ushort>()->default_value(100)->notifier(bind(&OptionsContainer::parseValue<ushort>, &options_container, "gibbs-burn-in", placeholders::_1)), "number of burn-in iterations.")
		("gibbs-samples", po::value<ushort>()->default_value(250)->notifier(bind(&OptionsContainer::parseValue<ushort>, &options_container, "gibbs-samples", placeholders::_1)), "number of Gibbs iterations.")
		("number-of-gibbs-chains", po::value<ushort>()->default_value(20)->notifier(bind(&OptionsContainer::parseValue<ushort>, &options_container, "number-of-gibbs-chains", placeholders::_1)), "number of parallel Gibbs sampling chains.")
		("kmer-subsampling-rate", po::value<float>()->default_value(0.1, "0.1")->notifier(bind(&OptionsContainer::parseValue<float>, &options_container, "kmer-subsampling-rate", placeholders::_1)), "subsampling rate for subsetting kmers used for genotype inference (a new subset is sampled for each Gibbs sampling chain).")
		("max-haplotype-variant-kmers", po::value<uint>()->default_value(500)->notifier(bind(&OptionsContainer::parseValue<uint>, &options_container, "max-haplotype-variant-kmers", placeholders::_1)), "maximum number of kmers used for genotype inference after subsampling across a haplotype candidate for each variant (a new subset is sampled for each Gibbs sampling chain).")
		("number-of-genomic-rate-gc-bias-bins", po::value<ushort>()->default_value(1)->notifier(bind(&OptionsContainer::parseValue<ushort>, &options_container, "number-of-genomic-rate-gc-bias-bins", placeholders::_1)), "number of genomic rate GC bias bins (a negative binomial distributed genomic rate is estimated for each bin).")
	;

	po::options_description noise("== Noise ==", 160);
	noise.add_options()

		("noise-rate-prior", po::value<string>()->default_value("1,1")->notifier(bind(&OptionsContainer::parseValuePair<double>, &options_container, "noise-rate-prior", placeholders::_1)), "gamma parameters for Poisson noise rate prior (<shape>,<scale>). All samples will use the same parameters.")
		("number-of-parameter-estimation-samples", po::value<ushort>()->default_value(100)->notifier(bind(&OptionsContainer::parseValue<ushort>, &options_container, "number-of-parameter-estimation-samples", placeholders::_1)), "number of parameter estimation iterations.")
		("number-of-parameter-estimation-SNVs", po::value<uint>()->default_value(1000000)->notifier(bind(&OptionsContainer::parseValue<uint>, &options_container, "number-of-parameter-estimation-SNVs", placeholders::_1)), "maximum number of autosomal SNVs to use for parameter estimation.")
	;

	po::options_description desc("## BayesTyper options ##");
	desc.add(help).add(required).add(general).add(graph).add(inference).add(noise);

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);

	if (vm.count("help") || argc == 1) {
		
		cout << desc << endl;
		return 1;
	}

	po::notify(vm);

	assert(Utils::ulong_overflow == (pow(2,64) - 1));
	assert(Utils::uint_overflow == (pow(2,32) - 1));
	assert(Utils::ushort_overflow == (pow(2,16) - 1));
	assert(Utils::uchar_overflow == (pow(2,8) - 1));

	assert((options_container.getValue<ushort>("kmer-size") == 31) or (options_container.getValue<ushort>("kmer-size") == 39) or (options_container.getValue<ushort>("kmer-size") == 47) or (options_container.getValue<ushort>("kmer-size") == 55) or (options_container.getValue<ushort>("kmer-size") == 63));
    assert(options_container.getValue<double>("copy-number-variant-threshold") >= 0);
    assert(options_container.getValue<double>("copy-number-variant-threshold") <= 1);
    assert(options_container.getValue<ushort>("max-number-of-sample-haplotype-candidates") > 0);
    assert(options_container.getValue<ushort>("max-number-of-sample-haplotype-candidates") < floor(Utils::ushort_overflow/float(30)));

    assert(options_container.getValue<ushort>("gibbs-burn-in") > 0);
    assert(options_container.getValue<ushort>("gibbs-samples") > 0);
    assert(options_container.getValue<ushort>("number-of-gibbs-chains") > 0);
    assert(options_container.getValue<float>("kmer-subsampling-rate") > 0);
    assert(options_container.getValue<float>("kmer-subsampling-rate") <= 1);
    assert(options_container.getValue<ushort>("number-of-genomic-rate-gc-bias-bins") <= options_container.getValue<ushort>("kmer-size"));

    assert((options_container.getValue<pair<double,double> >("noise-rate-prior").first) > 0);
    assert((options_container.getValue<pair<double,double> >("noise-rate-prior").second) > 0);
    assert(options_container.getValue<ushort>("number-of-parameter-estimation-samples") > 0);
    assert(options_container.getValue<uint>("number-of-parameter-estimation-SNVs") > 0);

	cout << "[" << Utils::getLocalTime() << "] Seeding pseudo-random number generator with " << options_container.getValue<int>("random-seed") << " ..." << endl;
	cout << "[" << Utils::getLocalTime() << "] Setting the kmer size to " << options_container.getValue<ushort>("kmer-size") << " ..." << endl;

	vector<Sample> samples;

	ifstream sample_file(options_container.getValue<string>("samples-file").c_str());
	assert(sample_file.is_open());

	string sample_line;

	while (sample_file.good()) {

		getline(sample_file, sample_line);

		if (sample_line.size() == 0) {

			continue;
		}

		vector<string> split_sample_line;
    	boost::split(split_sample_line, sample_line, boost::is_any_of("\t"));
    	assert(split_sample_line.size() == 3);

    	Sample sample;
    	sample.name = split_sample_line.at(0);

    	if (split_sample_line.at(1) == "M") {

    		sample.gender = Utils::Gender::Male;

    	} else {

    		assert(split_sample_line.at(1) == "F");
    		sample.gender = Utils::Gender::Female;
    	}

    	sample.file = split_sample_line.at(2);

		ifstream kmc_database_prefix_file(split_sample_line.at(2) + ".kmc_pre");
		ifstream kmc_database_suffix_file(split_sample_line.at(2) + ".kmc_suf");

		samples.push_back(sample);

		if (!kmc_database_prefix_file.good() or !kmc_database_suffix_file.good()) {

			cout << "\nERROR: " << split_sample_line.at(2) << ".kmc_pre or " << split_sample_line.at(2) << ".kmc_suf does not exist - or you do not have sufficient permissions.\n" << endl;
			exit(1);
		}
	}

    sample_file.close();

    if (samples.empty()) {

		cout << "\nERROR: Samples file empty.\n" << endl;
		exit(1);    	
    }

    if (samples.size() > 30) {

		cout << "\nERROR: The maximum number of samples supported by BayesTyper is currently 30.\n" << endl;
		exit(1);    	
    }

    const ushort num_samples = samples.size();

	cout << "\n[" << Utils::getLocalTime() << "] Parsed information for " << num_samples << " sample(s)\n" << endl;

	vector<VariantClusterGroup* > variant_cluster_groups; 

	KmerFactory kmer_factory(options_container);
	KmerHash * kmer_hash;

	switch (uchar(options_container.getValue<ushort>("kmer-size"))) {

		case 31: {

			kmer_hash = kmer_factory.initKmerHash<31>(&variant_cluster_groups, samples);
		    break;
		}

		case 39: {

			kmer_hash = kmer_factory.initKmerHash<39>(&variant_cluster_groups, samples);
		    break;
		}

		case 47: {

			kmer_hash = kmer_factory.initKmerHash<47>(&variant_cluster_groups, samples);
		    break;
		}

		case 55: {

			kmer_hash = kmer_factory.initKmerHash<55>(&variant_cluster_groups, samples);
		    break;
		}

		case 63: {

			kmer_hash = kmer_factory.initKmerHash<63>(&variant_cluster_groups, samples);
		    break;
		}

		default: {

			cout << "\nERROR: Invalid kmer size.\n" << endl;
			exit(1);
		}
	}

	InferenceEngine inference_engine(samples, options_container);
  
	CountDistribution count_distribution(num_samples, kmer_factory.genomicCountDistributions(), options_container);

	cout << "\n[" << Utils::getLocalTime() << "] " << Utils::getMaxMemoryUsage() << endl;

	inference_engine.estimateNoiseParameters(&count_distribution, &variant_cluster_groups, kmer_hash, samples);
	count_distribution.writeNoiseParameterEstimates(options_container.getValue<string>("output-prefix"), samples);

	cout << "\n[" << Utils::getLocalTime() << "] " << Utils::getMaxMemoryUsage() << endl;

    GenotypeWriter genotype_writer(samples, options_container.getValue<string>("output-prefix"), options_container.getValue<ushort>("threads"));
	
	inference_engine.genotypeVariantClusterGroups(&variant_cluster_groups, kmer_factory.numberOfVariants(), kmer_hash, count_distribution, samples, &genotype_writer);

	genotype_writer.completedGenotyping();

    genotype_writer.writeSampleAlleleKmerFractionCumDistFunc();
    genotype_writer.writeSampleAlleleKmerCoverageCumDistFunc();

	cout << "\n[" << Utils::getLocalTime() << "] " << Utils::getMaxMemoryUsage() << endl;

	delete kmer_hash;
    
	uint num_written_variants = genotype_writer.writeGenotypesToVariantCallFormat(options_container.getValue<string>("vcf-file"), options_container.getValue<Regions>("chromosome-regions"), options_container.writeHeader(), options_container.getValue<string>("genome-file"), kmer_factory.numberOfVariants());

	cout << "[" << Utils::getLocalTime() << "] Wrote " << num_written_variants << " variants" << endl;

	cout << "\n[" << Utils::getLocalTime() << "] BayesTyper completed succesfully!\n" << endl;

	return 0;
}
