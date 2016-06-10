
/*
main.cpp - This file is part of BayesTyper (v0.9)


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

	const string program_version = "v0.9";

	OptionsContainer options_container(program_version, Utils::getLocalTime());

	cout << "\n[" << Utils::getLocalTime() << "] " << "You are using BayesTyper (" << program_version << ")\n" << endl;
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

		("decoy-file,d", po::value<string>()->default_value("")->notifier(bind(&OptionsContainer::parseValue<string>, &options_container, "decoy-file", placeholders::_1)), "decoy sequences.")
		("chromosome-regions,c", po::value<string>()->default_value("")->notifier(bind(&OptionsContainer::parseRegions, &options_container, "chromosome-regions", placeholders::_1)), "infer genotypes only in the supplied list of chromosome regions (<chr1>:<start1>-<end1>,<chr2>,...); for whole chromosomes the start and end positions can be omitted. If no list is given the whole genome will be genotyped.")
    	("random-seed,r", po::value<uint>()->default_value(time(nullptr), "unix time")->notifier(bind(&OptionsContainer::parseValue<uint>, &options_container, "random-seed", placeholders::_1)), "seed for pseudo-random number generator (produces reproducible results only when running with one thread).")
		("kmer-size,k", po::value<ushort>()->default_value(default_kmer_size)->notifier(bind(&OptionsContainer::parseValue<ushort>, &options_container, "kmer-size", placeholders::_1)), "kmer size (31, 39, 47, 55 or 63).")
		("threads,p", po::value<ushort>()->default_value(1)->notifier(bind(&OptionsContainer::parseValue<ushort>, &options_container, "threads", placeholders::_1)), "number of threads used (+= 2 I/O threads).")
		("output-prefix,o", po::value<string>()->default_value("bayesTyper_genotypes")->notifier(bind(&OptionsContainer::parseValue<string>, &options_container, "output-prefix", placeholders::_1)), "prefix for output VCF file (suffix: \".vcf\").")
		("maximum-allele-length", po::value<uint>()->default_value(3000000)->notifier(bind(&OptionsContainer::parseValue<uint>, &options_container, "maximum-allele-length", placeholders::_1)), "exclude alleles (reference and alternative) longer than <value>.")
		("minimum-number-of-kmers-per-allele", po::value<ushort>()->default_value(0)->notifier(bind(&OptionsContainer::parseValue<ushort>, &options_container, "minimum-number-of-kmers-per-allele", placeholders::_1)), "exclude alleles with less than <value> kmers (<value> must be less than <kmer-size>).")
		("number-of-haplotype-candidates-per-sample", po::value<ushort>()->default_value(8)->notifier(bind(&OptionsContainer::parseValue<ushort>, &options_container, "number-of-haplotype-candidates-per-sample", placeholders::_1)), "maximum number of haplotype candidates per sample (dummy sample: <=24 candidates).")
	;

	po::options_description sampling("== Sampling ==", 160);
	sampling.add_options()

		("gibbs-burn-in,b", po::value<ushort>()->default_value(100)->notifier(bind(&OptionsContainer::parseValue<ushort>, &options_container, "gibbs-burn-in", placeholders::_1)), "number of burn-in iterations.")
		("gibbs-samples,i", po::value<ushort>()->default_value(100)->notifier(bind(&OptionsContainer::parseValue<ushort>, &options_container, "gibbs-samples", placeholders::_1)), "number of gibbs iterations.")
		("number-of-gibbs-chains,n", po::value<ushort>()->default_value(10)->notifier(bind(&OptionsContainer::parseValue<ushort>, &options_container, "number-of-gibbs-chains", placeholders::_1)), "number of parallel gibbs sampling chains.")
		("number-of-parameter-estimation-samples,e", po::value<ushort>()->default_value(100)->notifier(bind(&OptionsContainer::parseValue<ushort>, &options_container, "number-of-parameter-estimation-samples", placeholders::_1)), "number of parameter estimation iterations.")
		("number-of-parameter-estimation-variants", po::value<uint>()->default_value(1000000)->notifier(bind(&OptionsContainer::parseValue<uint>, &options_container, "number-of-parameter-estimation-variants", placeholders::_1)), "expected number of variants to use for parameter estimation.")
		("maximum-number-of-multicluster-kmers", po::value<uint>()->default_value(10000)->notifier(bind(&OptionsContainer::parseValue<uint>, &options_container, "maximum-number-of-multicluster-kmers", placeholders::_1)), "maximum number of multicluster kmers per variant cluster used when gibbs sampling (a new subset is sampled for each chain equally among all haplotype candidates).")
	;

	po::options_description noise("== Noise ==", 160);
	noise.add_options()

		("noise-zero-inflation-priors", po::value<string>()->default_value("1,1")->notifier(bind(&OptionsContainer::parsePriors, &options_container, "noise-zero-inflation-priors", placeholders::_1)), "beta parameters for noise zero inflation prior(s) (<alpha1>,<beta1>:<alpha2>,<beta2>:...); a pair for each noise component. All samples will use the same set of priors.")
		("noise-rate-priors", po::value<string>()->default_value("1,1")->notifier(bind(&OptionsContainer::parsePriors, &options_container, "noise-rate-priors", placeholders::_1)), "gamma parameters for noise rate prior(s) (<shape1>,<scale1>:<shape2>,<scale2>:...); a pair for each noise component. All samples will use the same set of priors.")
	;

	po::options_description filter("== Filter ==", 160);
	filter.add_options()

		("minimum-number-of-observed-informative-kmers", po::value<ushort>()->default_value(1)->notifier(bind(&OptionsContainer::parseValue<ushort>, &options_container, "minimum-number-of-observed-informative-kmers", placeholders::_1)), "Filter sampled alleles with less than <value> observed (coverage above zero) informative (unique multiplicity compared to the other alleles in the variant) kmers. Missing \"*\" and not sampled alleles (values of -1) is not filtered.")
	;

	po::options_description desc("## BayesTyper options ##");
	desc.add(help).add(required).add(general).add(sampling).add(noise).add(filter);

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

	assert(options_container.getValue<ushort>("minimum-number-of-kmers-per-allele") < options_container.getValue<ushort>("kmer-size"));
    assert(options_container.getValue<ushort>("number-of-haplotype-candidates-per-sample") > 0);

    assert(options_container.getValue<ushort>("gibbs-burn-in") > 0);
    assert(options_container.getValue<ushort>("gibbs-samples") > 0);
    assert(options_container.getValue<ushort>("number-of-gibbs-chains") > 0);
    assert(options_container.getValue<ushort>("number-of-parameter-estimation-samples") > 0);
    assert(options_container.getValue<ushort>("number-of-parameter-estimation-variants") > 0);

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

    		sample.sex = Utils::Sex::Male;

    	} else {

    		assert(split_sample_line.at(1) == "F");
    		sample.sex = Utils::Sex::Female;
    	}

    	sample.file = split_sample_line.at(2);

		ifstream kmc_database_prefix_file(split_sample_line.at(2) + ".kmc_pre");
		ifstream kmc_database_suffix_file(split_sample_line.at(2) + ".kmc_suf");

		samples.push_back(sample);

		if (!kmc_database_prefix_file.good() or !kmc_database_suffix_file.good()) {

			cout << "\n ERROR: " << split_sample_line.at(2) << ".kmc_pre or " << split_sample_line.at(2) << ".kmc_suf does not exist - or you do not have sufficient permissions." << endl;
			exit(1);
		}
	}

    sample_file.close();

    assert(!(samples.empty()));
    assert(samples.size() < Utils::ushort_overflow);

    const ushort num_samples = samples.size();

    assert(static_cast<uint>(options_container.getValue<ushort>("number-of-haplotype-candidates-per-sample") * num_samples) <= static_cast<uint>(Utils::ushort_overflow - 2));

    assert(static_cast<uint>(options_container.getValue<vector<pair<double,double> > >("noise-zero-inflation-priors").size()) < Utils::uchar_overflow);

	const uchar num_noise_sources = options_container.getValue<vector<pair<double,double> > >("noise-zero-inflation-priors").size();

	assert(num_noise_sources == static_cast<uchar>(options_container.getValue<vector<pair<double,double> > >("noise-rate-priors").size()));
	assert(num_noise_sources <= 8);

	cout << "\n[" << Utils::getLocalTime() << "] Parsed parameters for " << num_samples << " sample(s) and " << to_string(num_noise_sources) << " noise source(s)\n" << endl;

	vector<VariantClusterGroup* > variant_cluster_groups; 

	KmerFactory kmer_factory(options_container);
	KmerHash * kmer_hash;

	switch (uchar(options_container.getValue<ushort>("kmer-size"))) {

		default: {

			cout << "\nERROR: Invalid kmer size.\n" << endl;

			exit(1);
			break;
		}

		case 31: {

			kmer_hash = kmer_factory.buildKmerHash<31>(&variant_cluster_groups, samples);
		    break;
		}

		case 39: {

			kmer_hash = kmer_factory.buildKmerHash<39>(&variant_cluster_groups, samples);
		    break;
		}

		case 47: {

			kmer_hash = kmer_factory.buildKmerHash<47>(&variant_cluster_groups, samples);
		    break;
		}

		case 55: {

			kmer_hash = kmer_factory.buildKmerHash<55>(&variant_cluster_groups, samples);
		    break;
		}

		case 63: {

			kmer_hash = kmer_factory.buildKmerHash<63>(&variant_cluster_groups, samples);
		    break;
		}
	}

    GenotypeWriter genotype_writer(samples, options_container);
	InferenceEngine inference_engine(samples, num_noise_sources, options_container);
  
	auto count_distribution = CountDistribution(num_samples, num_noise_sources, kmer_factory.genomicCountDistributions(), options_container);

	cout << "\n[" << Utils::getLocalTime() << "] " << Utils::getMaxMemoryUsage() << endl;

	cout << "\n[" << Utils::getLocalTime() << "] Estimating noise model parameters ..." << endl;

	inference_engine.estimateNoiseParameters(&count_distribution, &variant_cluster_groups, kmer_hash, samples);
	count_distribution.writeParameterSamples(samples);

	cout << "\n[" << Utils::getLocalTime() << "] Caching model probabilities ..." << endl;
    
	count_distribution.updateMaxMultiplicityCache(Utils::uchar_overflow);

    cout << "[" << Utils::getLocalTime() << "] Finished caching\n" << endl;

    cout << "[" << Utils::getLocalTime() << "] Sorting variant cluster groups by decreasing complexity (number of kmers and haplotype candidates) ..." << endl;

    sort(variant_cluster_groups.begin(), variant_cluster_groups.end(), VariantClusterGroupCompare);

    cout << "[" << Utils::getLocalTime() << "] Finished sorting\n" << endl;

	cout << "\n[" << Utils::getLocalTime() << "] " << Utils::getMaxMemoryUsage() << endl;

	cout << "\n[" << Utils::getLocalTime() << "] Genotyping variant cluster groups ..." << endl;
	
	inference_engine.genotypeVariantClusterGroups(&variant_cluster_groups, kmer_hash, count_distribution, samples, &genotype_writer);

	genotype_writer.completedGenotyping();

	cout << "\n[" << Utils::getLocalTime() << "] " << Utils::getMaxMemoryUsage() << endl;

	delete kmer_hash;

	cout << "\n[" << Utils::getLocalTime() << "] Writing genotypes to variant file " << options_container.getValue<string>("output-prefix") << ".vcf ..." << endl;
    
	uint num_written_variants = genotype_writer.writeGenotypesToVariantFile(options_container.getValue<string>("vcf-file"), options_container.getValue<Regions>("chromosome-regions"), samples, options_container.writeVCFFileHeader(), options_container.getValue<string>("genome-file"), kmer_factory.numberOfVariants());

	cout << "[" << Utils::getLocalTime() << "] Wrote " << num_written_variants << " variants (" << kmer_factory.numberOfVariants() - num_written_variants << " were unsupported or skipped)" << endl;

	cout << "\n[" << Utils::getLocalTime() << "] BayesTyper completed succesfully!!!\n" << endl;

	return 0;
}
