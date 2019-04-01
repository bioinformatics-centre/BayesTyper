
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


#include <vector>
#include <string>
#include <iostream>
#include <sstream>

#include "boost/program_options.hpp"

#include "Utils.hpp"
#include "CompareOperators.hpp"

#include "MakeBloom.hpp"
#include "ConvertAllele.hpp"
#include "Combine.hpp"
#include "Filter.hpp"
#include "Annotate.hpp"
#include "AddAttributes.hpp"

namespace po = boost::program_options;

// static const string parents_trio_regex_default = "^[0-9]+((-01)|(-02)){1}$";

using namespace std;

string generateVariantFileOutput(const string & output_prefix, const bool is_compressed) {

	stringstream suffix;
	suffix << output_prefix << ".vcf";

	if (is_compressed) {

		suffix << ".gz";
	}

	return suffix.str();
}

int main (int argc, char * const argv[]) {

	std::cout << "\n[" << Utils::getLocalTime() << "] " << "You are using BayesTyperTools (" << BT_VERSION << ")\n" << std::endl;

	stringstream command_info;

	command_info << "Usage: bayesTyperTools <command> [options]" << endl;
	command_info << "\nCommands:\n" << endl;
	command_info << "\tmakeBloom\t\tcreate kmer bloom filter" << endl;
	command_info << "\tconvertAllele\t\tconvert allele IDs to sequence" << endl;
	command_info << "\tcombine\t\t\tcombine callsets (vertical)" << endl;
	command_info << "\tfilter\t\t\tfilter variants, alleles and/or samples" << endl;
	command_info << "\tannotate\t\tannotate alleles" << endl;
	command_info << "\taddAttributes\t\tadd variant, allele and/or trio attributes" << endl;

	if (argc == 1) {

		cout << command_info.str() << endl;

	} else {

		po::options_description help_options("", 160);
		help_options.add_options()

		   ("help,h", "produce help message for options")
		;

		bool gzip_output;

		po::options_description general_options("== General ==", 160);
		general_options.add_options()

			("gzip-output,z", po::value<bool>(&gzip_output)->default_value(false)->implicit_value(true), "compress output file(s) using gzip.")
		;

		if (strcmp(argv[1],"makeBloom") == 0) {

			string kmc_table_prefix;
			
			ushort num_threads;
			bool run_test = false;

			float false_positive_rate;

			po::options_description required_options("== Required ==", 160);
			required_options.add_options()

				("kmc-table-prefix,k", po::value<string>(&kmc_table_prefix)->required(), "KMC kmer table prefix. Output is written as <kmc-table-prefix>.bloomMeta and <kmc-table-prefix>.bloomData.")
			;

			po::options_description general_options_makebloom("== General ==", 160);
			general_options_makebloom.add_options()

				("num-threads,p", po::value<ushort>(&num_threads)->default_value(1), "number of threads used (+= 1 I/O thread).")

				// ("run-test", po::value<bool>(&run_test)->default_value(false)->implicit_value(true), "test bloom filter. WARNING: Memory intensive!")
			;

			po::options_description parameters_options("== Parameters ==", 160);
			parameters_options.add_options()

				("false-positive-rate", po::value<float>(&false_positive_rate)->default_value(0.001, "0.001"), "bloom filter false positive rate.")

			;

			po::options_description desc("## BayesTyperTools makeBloom ##");
			desc.add(help_options).add(required_options).add(general_options_makebloom).add(parameters_options);

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);

			if (vm.count("help") || argc == 2) {
			   cout << desc << endl;
			   return 1;
			}

			po::notify(vm);

			MakeBloom make_bloom;
			make_bloom.kmc2bloomFile(kmc_table_prefix, false_positive_rate, run_test, num_threads);

		} else if (strcmp(argv[1],"convertAllele") == 0) {

			string variant_file;
			string genome_file;
			string output_prefix;

			string mei_file;
			string alt_file;
			bool keep_imprecise;
			bool keep_partial;

			po::options_description required_options("== Required ==", 160);
			required_options.add_options()

				("variant-file,v", po::value<string>(&variant_file)->required(), "variant file (vcf format).")
				("genome-file,g", po::value<string>(&genome_file)->required(), "reference genome file (fasta format).")
				("output-prefix,o", po::value<string>(&output_prefix)->required(), "output prefix.")
			;

			po::options_description alleles_options("== Alleles ==", 160);
			alleles_options.add_options()

				("alt-file", po::value<string>(&alt_file)->default_value(""), "alternative allele file (fasta format). Sequence name in fasta (>\"name\") should match <\"name\">.")
				("mei-file", po::value<string>(&mei_file)->default_value(""), "mobile element insertion(s) file (fasta format). Sequence name in fasta (>\"name\") should match <INS:ME:\"name\">.")
				("keep-imprecise", po::value<bool>(&keep_imprecise)->default_value(false)->implicit_value(true), "do not filter imprecise variants")
				("keep-partial", po::value<bool>(&keep_partial)->default_value(false)->implicit_value(true), "keep partial insertions where the center and length is unknown (Manta output supported). The known left and right side is connected with ten N's.")
			;

			po::options_description desc("## BayesTyperTools convertAllele ##");
			desc.add(help_options).add(required_options).add(general_options).add(alleles_options);

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);

			if (vm.count("help") || argc == 2) {
			   cout << desc << endl;
			   return 1;
			}

			po::notify(vm);

			ConvertAllele::convertAllele(variant_file, genome_file, generateVariantFileOutput(output_prefix, gzip_output), alt_file, mei_file, keep_imprecise, keep_partial);

		} else if (strcmp(argv[1],"combine") == 0) {

			string variant_files;
			string output_prefix;

			bool filter_ambiguous_alleles;

			po::options_description required_options("== Required ==", 160);
			required_options.add_options()

				("variant-files,v", po::value<string>(&variant_files)->required(), "comma-separated list of name and variant file (vcf format) pairs (<name>:<file>).")
				("output-prefix,o", po::value<string>(&output_prefix)->required(), "output prefix.")
			;

			po::options_description filters_options("== Filters ==", 160);
			filters_options.add_options()

				("filter-ambiguous-alleles", po::value<bool>(&filter_ambiguous_alleles)->default_value(false)->implicit_value(true), "filter alleles (including reference) with ambiguous nucleotides (non ACGT).")
			;

			po::options_description desc("## BayesTyperTools combine ##");
			desc.add(help_options).add(required_options).add(general_options).add(filters_options);

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);

			if (vm.count("help") || argc == 2) {
			   cout << desc << endl;
			   return 1;
			}

			po::notify(vm);

			Combine::combine(Utils::splitString(variant_files, ','), generateVariantFileOutput(output_prefix, gzip_output), filter_ambiguous_alleles);
		
		} else if (strcmp(argv[1],"filter") == 0) {

			string variant_file;
			string output_prefix;

			string kmer_coverage_file;
	        Filter::FilterValues filter_values;

			po::options_description required_options("== Required ==", 160);
			required_options.add_options()

				("variant-file,v", po::value<string>(&variant_file)->required(), "variant file (vcf format).")
				("output-prefix,o", po::value<string>(&output_prefix)->required(), "output prefix.")
			;

			po::options_description filters_options("== Filters ==", 160);
			filters_options.add_options()

				("min-homozygote-genotypes", po::value<uint>(&filter_values.min_homozygote)->default_value(0), "filter variants with less than <value> homozygote genotypes (calculated before other filters).")
				("min-genotype-posterior", po::value<float>(&filter_values.min_gpp_value)->default_value(0.99, "0.99"), "filter genotypes with a posterior probability (GPP) below <value>.")
				("min-number-of-kmers", po::value<float>(&filter_values.min_nak_value)->default_value(1, "1"), "filter sampled alleles with less than <value> kmers (NAK).")
				("kmer-coverage-file", po::value<string>(&kmer_coverage_file)->default_value("bayestyper_genomic_parameters.txt"), "sample kmer coverage file used for filtering sampled alleles with a low fraction of observed kmers (FAK).")
			;

			po::options_description desc("## BayesTyperTools filter ##");
			desc.add(help_options).add(required_options).add(general_options).add(filters_options);

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);

			if (vm.count("help") || argc == 2) {
			   cout << desc << endl;
			   return 1;
			}

			po::notify(vm);

			Filter::filter(variant_file, generateVariantFileOutput(output_prefix, gzip_output), kmer_coverage_file, filter_values);

		} else if (strcmp(argv[1],"annotate") == 0) {

			string variant_file;
			string annotation_file;
			string output_prefix;

			bool clear_prev_annotation;

			float match_threshold;
			float window_size_scale;

			po::options_description required_options("== Required ==", 160);
			required_options.add_options()

				("variant-file,v", po::value<string>(&variant_file)->required(), "variant file (vcf format).")
				("annotation-file,a", po::value<string>(&annotation_file)->required(), "annotation file (vcf format).")
				("output-prefix,o", po::value<string>(&output_prefix)->required(), "output prefix.")
			;

			auto general_options_annotate = general_options;
			general_options_annotate.add_options()

				("clear-prev-annotation,c", po::value<bool>(&clear_prev_annotation)->default_value(false)->implicit_value(true), "clear previous annotations (variant id and AAI).")
			;

			po::options_description parameters_options("== Parameters ==", 160);
			parameters_options.add_options()

				("match-threshold", po::value<float>(&match_threshold)->default_value(0.5, "0.5"), "minimum sequence overlap between input allele and annotation allele.")
				("window-size-scale", po::value<float>(&window_size_scale)->default_value(3, "3"), "window size allele length scaling factor.")
			;

			po::options_description desc("## BayesTyperTools annotate ##");
			desc.add(help_options).add(required_options).add(general_options_annotate).add(parameters_options);

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);

			if (vm.count("help") || argc == 2) {

			   	cout << desc << endl;
				return 1;
			}

			po::notify(vm);

			Annotate::annotate(variant_file, annotation_file, generateVariantFileOutput(output_prefix, gzip_output), match_threshold, window_size_scale, clear_prev_annotation);

		} else if (strcmp(argv[1],"addAttributes") == 0) {

			string variant_file;
			string output_prefix;

			string genome_file;
			string repeat_file;
			string indepedent_samples_regex_str;
			string trio_sample_info_str;

			po::options_description required_options("== Required ==", 160);
			required_options.add_options()

				("variant-file,v", po::value<string>(&variant_file)->required(), "variant file (vcf format).")
				("output-prefix,o", po::value<string>(&output_prefix)->required(), "output prefix.")
			;

			po::options_description attributes_options("== Attributes ==", 160);
			attributes_options.add_options()

				("genome-file", po::value<string>(&genome_file), "reference genome file (fasta format) used for homopolymer length (HPL) calculation. If not specified HPL will not be calculated.")
				("repeat-file", po::value<string>(&repeat_file), "repeatmasker file used for repeat annotation (RMA). If not specified RMA will not be annotated.")
				("independent-samples-regex", po::value<string>(&indepedent_samples_regex_str), "regular expression for matching independent samples (e.g. parents in a trio) used for absolute inbreeding coefficient (IBC) calculation. If not specified IBC will not be calculated.")
				("trio-sample-info", po::value<string>(&trio_sample_info_str), "trio sample id information used for concordance (CONC) calculation (<father>,<mother>,<child>:<father>,<mother>,<child>:...). If not specified CONC will not be calculated.")
			;

			po::options_description desc("## BayesTyperTools addAttributes ##");
			desc.add(help_options).add(required_options).add(general_options).add(attributes_options);

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);

			if (vm.count("help") || argc == 2) {
			   cout << desc << endl;
			   return 1;
			}

			po::notify(vm);

			AddAttributes::addAttributes(variant_file, generateVariantFileOutput(output_prefix, gzip_output), genome_file, repeat_file, indepedent_samples_regex_str, trio_sample_info_str);

		} else {

			cout << command_info.str() << endl;
		}
	}

	return 0;
}
