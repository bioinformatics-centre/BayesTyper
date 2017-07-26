
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


#include <vector>
#include <string>
#include <iostream>
#include <sstream>

#include "boost/program_options.hpp"

#include "Utils.hpp"
#include "CompareOperators.hpp"

#include "convertAlleleId.hpp"
#include "combine.hpp"
#include "split.hpp"
#include "merge.hpp"
#include "writeIndels.hpp"
#include "annotate.hpp"
#include "addAttributes.hpp"
#include "filter.hpp"
#include "getSummary.hpp"

namespace po = boost::program_options;

// static const string parents_trio_regex_default = "^[0-9]+((-01)|(-02)){1}$";
static const string parents_trio_regex_default = "^$";

using namespace std;

int main (int argc, char * const argv[]) {

	std::cout << "\n[" << Utils::getLocalTime() << "] " << "You are using BayesTyperTools (" << BT_VERSION << ")\n" << std::endl;

	stringstream command_info;

	command_info << "Usage BayesTyperTools <command> [options]" << endl;
	command_info << "\nCommands:\n" << endl;
	command_info << "\tconvertAlleleId\t\tconvert allele IDs to sequence" << endl;
	command_info << "\tcombine\t\t\tcombine callsets (vertical)" << endl;
	command_info << "\tsplit\t\t\tsplit callset(s)" << endl;
	command_info << "\tmerge\t\t\tmerge callsets (horizontal)" << endl;
	command_info << "\tannotate\t\tannotate alleles" << endl;
	command_info << "\twriteIndels\t\twrite indel sequences" << endl;
	command_info << "\taddAttributes\t\tadd variant, allele and/or trio attributes" << endl;
	command_info << "\tfilter\t\t\tfilter variants, alleles and/or samples" << endl;
	command_info << "\tgetSummary\t\tget variant, allele, trio and sample summary" << endl;

	if (argc == 1) {

		cout << command_info.str() << endl;

	} else {

		if (strcmp(argv[1],"convertAlleleId") == 0) {

			string vcf_filename;
			string genome_filename;
			string output_prefix;

			string mei_filename = "";
			bool keep_imprecise_variants;

			po::options_description general("", 150);
			general.add_options()
			   ("help,h", "produce help message for base options")
			;

			po::options_description required_options("== Required options ==", 160);
			required_options.add_options()

				("vcf-filename,v", po::value<std::string>(&vcf_filename)->required(), "variant file.")
				("genome-filename,g", po::value<std::string>(&genome_filename)->required(), "genome file.")
				("output-prefix,o", po::value<std::string>(&output_prefix)->required(), "output prefix (suffix: \".vcf\").")
			;

			po::options_description optional_options("== Optional options ==", 160);
			optional_options.add_options()

				("mobile-element-insertion-filename", po::value<std::string>(&mei_filename), "fasta file contaning mobile element insertion sequence(s) (\"header\" should be a perfect match to the <INS:ME:\"header\">)")
				("keep-imprecise-variants", po::value<bool>(&keep_imprecise_variants)->default_value(false)->implicit_value(true), "do not filter imprecise variants")
			;

			po::options_description desc("## BayesTyperTools convertAlleleId ##");
			desc.add(general).add(required_options).add(optional_options);

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);

			if (vm.count("help") || argc == 2) {
			   cout << desc << endl;
			   return 1;
			}

			po::notify(vm);

			ConvertAlleleId::convertAlleleId(vcf_filename, genome_filename, output_prefix, mei_filename, keep_imprecise_variants);

		} else if (strcmp(argv[1],"combine") == 0) {

			string vcf_filenames_str;
			string output_prefix;

			bool exclude_ambiguous_alleles;

			po::options_description general("", 150);
			general.add_options()
			   ("help,h", "produce help message for base options")
			;

			po::options_description required_options("== Required options ==", 160);
			required_options.add_options()

				("vcf-filenames,v", po::value<std::string>(&vcf_filenames_str)->required(), "call-set names and file (<call-set>:<filename>,<call-set>:<filename>,...).")
				("output-prefix,o", po::value<std::string>(&output_prefix)->required(), "output prefix (suffix: \".vcf\").")
			;

			po::options_description optional_options("== Optional options ==", 160);
			optional_options.add_options()

				("exclude-ambiguous-alleles", po::value<bool>(&exclude_ambiguous_alleles)->default_value(false)->implicit_value(true), "exclude alleles (including reference) with ambiguous nucleotides (non ACGT).")
			;

			po::options_description desc("## BayesTyperTools combine ##");
			desc.add(general).add(required_options).add(optional_options);

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);

			if (vm.count("help") || argc == 2) {
			   cout << desc << endl;
			   return 1;
			}

			po::notify(vm);

			Combine::combine(Utils::splitString(vcf_filenames_str, ','), output_prefix, exclude_ambiguous_alleles);

		} else if (strcmp(argv[1],"split") == 0) {

			string vcf_filenames_str;
			uint min_batch_size;
			bool use_vcg_split;
			bool use_chr_split;

			po::options_description general("", 150);
			general.add_options()
			   ("help,h", "produce help message for base options")
			;

			po::options_description required_options("== Required options ==", 160);
			required_options.add_options()

				("vcf-filenames,v", po::value<std::string>(&vcf_filenames_str)->required(), "batch name and file (<batch>:<filename>,<batch>:<filename>,...). The variants need to be identical between input batches. The number of threads will be equal to the number of input batches.")
			;

			po::options_description optional_options("== Optional options ==", 160);
			optional_options.add_options()

				("min-batch-size", po::value<uint>(&min_batch_size)->default_value(50000), "minimum number of variants in a variant batch.")
				("use-vcg-split", po::value<bool>(&use_vcg_split)->default_value(false)->implicit_value(true), "require split between variant cluster groups.")
				("use-chr-split", po::value<bool>(&use_chr_split)->default_value(false)->implicit_value(true), "require split between chromosomes (overwrites <min-batch-size>).")
			;

			po::options_description desc("## BayesTyperTools split ##");
			desc.add(general).add(required_options).add(optional_options);

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);

			if (vm.count("help") || argc == 2) {
			   cout << desc << endl;
			   return 1;
			}

			po::notify(vm);

			assert(!use_vcg_split or !use_chr_split);

			if (use_chr_split) {

				min_batch_size = 0;
			}

			Split::split(Utils::splitString(vcf_filenames_str, ','), min_batch_size, use_vcg_split, use_chr_split);

		} else if (strcmp(argv[1],"merge") == 0) {

			string vcf_filenames_str;
			string output_prefix;
			string parents_trio_regex;

			po::options_description general("", 150);
			general.add_options()
			   ("help,h", "produce help message for base options")
			;

			po::options_description required_options("== Required options ==", 160);
			required_options.add_options()

				("vcf-filenames,v", po::value<std::string>(&vcf_filenames_str)->required(), "variant files.")
				("output-prefix,o", po::value<std::string>(&output_prefix)->required(), "output prefix (suffix: \".vcf\").")
			;

			po::options_description desc("## BayesTyperTools merge ##");
			desc.add(general).add(required_options);

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);

			if (vm.count("help") || argc == 2) {
			   cout << desc << endl;
			   return 1;
			}

			po::notify(vm);

			Merge::merge(Utils::splitString(vcf_filenames_str, ','), output_prefix);

		} else if (strcmp(argv[1],"annotate") == 0) {

			string vcf_filename;
			string annotation_filename;
			string output_prefix;
			
			float match_threshold;
			float window_size_scale;
			bool overwrite_prev_anno;

			po::options_description general("", 150);
			general.add_options()
			   ("help,h", "produce help message for base options")
			;

			po::options_description required_options("== Required options ==", 160);
			required_options.add_options()

				("vcf-filename,v", po::value<std::string>(&vcf_filename)->required(), "variant file.")
				("annotation,a", po::value<std::string>(&annotation_filename)->required(), "annotation file (vcf format). Needs to be sorted with contig order in the header identical to <vcf-filename>.")
				("output-prefix,o", po::value<std::string>(&output_prefix)->required(), "output prefix (suffix: \".vcf\").")
			;

			po::options_description optional_options("== Optional options ==", 160);
			optional_options.add_options()

				("match-threshold", po::value<float>(&match_threshold)->default_value(0.5, "0.5"), "minimum sequence overlap between input allele and annotation allele.")
				("window-size-scale", po::value<float>(&window_size_scale)->default_value(3, "3"), "window size allele length scaling factor.")
				("overwrite-prev-annotation", po::value<bool>(&overwrite_prev_anno)->default_value(false)->implicit_value(true), "overwrite previous annotations (variant id and \"AAI\" attribute).")
			;

			po::options_description desc("## BayesTyperTools annotate ##");
			desc.add(general).add(required_options).add(optional_options);

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);

			if (vm.count("help") || argc == 2) {

			   	cout << desc << endl;
				return 1;
			}

			po::notify(vm);

			Annotate::annotate(vcf_filename, annotation_filename, output_prefix, match_threshold, window_size_scale, overwrite_prev_anno);

		} else if (strcmp(argv[1],"writeIndels") == 0) {

			string vcf_filename;
			string output_prefix;
			uint min_indel_length;

			po::options_description general("", 150);
			general.add_options()
			   ("help,h", "produce help message for base options")
			;

			po::options_description required_options("== Required options ==", 160);
			required_options.add_options()

				("vcf-filename,v", po::value<std::string>(&vcf_filename)->required(), "variant file.")
				("output-prefix,o", po::value<std::string>(&output_prefix)->required(), "output prefix (suffix: \"_indelseqs.fa\")")
			;

			po::options_description optional_options("== Optional options ==", 160);
			optional_options.add_options()

				("min-indel-length", po::value<uint>(&min_indel_length)->default_value(10), "minimum indel length")
			;

			po::options_description desc("## BayesTyperTools writeIndels ##");
			desc.add(general).add(required_options).add(optional_options);

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);

			if (vm.count("help") || argc == 2) {
			   cout << desc << endl;
			   return 1;
			}

			po::notify(vm);
			WriteIndels::writeIndels(vcf_filename, output_prefix, min_indel_length);

		} else if (strcmp(argv[1],"addAttributes") == 0) {

			string vcf_filename;
			string output_prefix;

			string genome_filename;
			string repeat_filename;
			string parents_trio_regex;
			string trio_info_str;
			string fak_cdf_filename;
			string mac_cdf_filename;

			po::options_description general("", 150);
			general.add_options()
			   ("help,h", "produce help message for base options")
			;

			po::options_description required_options("== Required options ==", 160);
			required_options.add_options()

				("vcf-filename,v", po::value<std::string>(&vcf_filename)->required(), "variant file.")
				("output-prefix,o", po::value<std::string>(&output_prefix)->required(), "output prefix (suffix: \".vcf\").")
			;

			po::options_description optional_options("== Optional options ==", 160);
			optional_options.add_options()

				("genome-filename", po::value<std::string>(&genome_filename), "genome file used for homopolymer length (HPL) calculation. If not specified HPL will not be calculated.")
				("repeatmasker-filename", po::value<string>(&repeat_filename), "repeatmasker file used for repeat annotation (RMA). If not specified RMA will not be annotated.")
				("parents-trio-reg-expression", po::value<string>(&parents_trio_regex)->default_value(parents_trio_regex_default), "regular expression for matching parents in trio samples used for absolute inbreeding coefficient (IBC) calculation. If not specified IBC will not be calculated.")
				("trio-info", po::value<std::string>(&trio_info_str), "trio sample id information used for concordance (CONC) calculation (<father>,<mother>,<child>:<father>,<mother>,<child>:...). If not specified CONC will not be calculated.")
				("fak-cdf-filename", po::value<std::string>(&fak_cdf_filename), "cumulative distribution function file for fraction of observed kmers quantile (FAKQ) calculation. If not specified FAKQ will not be calculated.")
				("mac-cdf-filename", po::value<std::string>(&mac_cdf_filename), "cumulative distribution function file for kmer coverage quantile (MACQ) calculation. If not specified MACQ will not be calculated.")
			;

			po::options_description desc("## BayesTyperTools addAttributes ##");
			desc.add(general).add(required_options).add(optional_options);

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);

			if (vm.count("help") || argc == 2) {
			   cout << desc << endl;
			   return 1;
			}

			po::notify(vm);

			AddAttributes::addAttributes(vcf_filename, output_prefix, genome_filename, repeat_filename, parents_trio_regex, trio_info_str, fak_cdf_filename, mac_cdf_filename);

		} else if (strcmp(argv[1],"filter") == 0) {

			string vcf_filename;
			string genome_filename;
			string output_prefix;
			string kmer_coverage_filename;

			string parents_trio_regex;
	        Filter::FilterValues filter_values;

			po::options_description general("", 150);
			general.add_options()
			   ("help,h", "produce help message for base options")
			;

			po::options_description required_options("== Required options ==", 160);
			required_options.add_options()

				("vcf-filename,v", po::value<std::string>(&vcf_filename)->required(), "variant file.")
				("genome-filename,g", po::value<std::string>(&genome_filename)->required(), "genome file.")
				("output-prefix,o", po::value<std::string>(&output_prefix)->required(), "output prefix (suffix: \".vcf\").")
			;

			po::options_description optional_options("== Optional options ==", 160);
			optional_options.add_options()

				("parents-trio-reg-expression", po::value<string>(&parents_trio_regex)->default_value(parents_trio_regex_default), "regular expression for matching parents in trio samples.")
				("max-inbreeding-coef", po::value<float>(&filter_values.max_inbreeding_coef)->default_value(1, "1"), "filter variants with an absolute inbreeding coefficient above <value> (calculated before filtering). Minimum 10 independent samples required for this filter.")
				("max-homopolymer-length", po::value<float>(&filter_values.max_homopolymer_length)->default_value(9, "9"), "filter insertions and deletions located in a homopolymer longer than <value>.")
				("min-number-of-kmers", po::value<float>(&filter_values.min_nak_value)->default_value(10, "10"), "filter sampled alleles with less than <value> kmers.")
				("kmer-coverage-filename", po::value<string>(&kmer_coverage_filename), "sample kmer coverage file used for filtering sampled alleles with a low fraction of observed kmers.")
				("min-allelic-balance", po::value<float>(&filter_values.min_allelic_balance)->default_value(0.3, "0.3"), "filter genotypes with an allelic balance below <value>.")
				("min-genotype-posterior", po::value<float>(&filter_values.min_gpp_value)->default_value(0.99, "0.99"), "filter genotypes with a posterior probability below <value>.")
				("min-number-of-homozygote", po::value<uint>(&filter_values.min_homozygote)->default_value(1), "filter variants with less than <value> homozygote unfiltered genotypes (zeroploid and haploid genotypes are counted as homozygote). Minimum 10 samples required for this filter.")
			;

			po::options_description desc("## BayesTyperTools filter ##");
			desc.add(general).add(required_options).add(optional_options);

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);

			if (vm.count("help") || argc == 2) {
			   cout << desc << endl;
			   return 1;
			}

			po::notify(vm);

			Filter::filter(vcf_filename, genome_filename, output_prefix, kmer_coverage_filename, parents_trio_regex, filter_values);

		} else if (strcmp(argv[1],"getSummary") == 0) {

			string vcf_filename;
			string output_prefix;

			string trio_info_str;

			po::options_description general("", 150);
			general.add_options()
			   ("help,h", "produce help message for base options")
			;

			po::options_description required_options("== Required options ==", 160);
			required_options.add_options()

				("vcf-filename,v", po::value<std::string>(&vcf_filename)->required(), "variant file.")
				// ("output-prefix,o", po::value<std::string>(&output_prefix)->required(), "output prefix (suffix: \"_variant.txt\", \"_allele.txt\", \"_trio.txt\" & \"_sample.txt\").")
				("output-prefix,o", po::value<std::string>(&output_prefix)->required(), "output prefix (suffix: \"_variant.txt\" & \"_allele.txt\"")
			;

			po::options_description optional_options("== Optional options ==", 160);
			optional_options.add_options()

				("trio-info", po::value<std::string>(&trio_info_str), "trio sample id information used for trio summary output (<father>,<mother>,<child>:<father>,<mother>,<child>:...).")
			;

			po::options_description desc("## BayesTyperTools getSummary ##");
			desc.add(general).add(required_options).add(optional_options);

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);

			if (vm.count("help") || argc == 2) {
			   cout << desc << endl;
			   return 1;
			}

			po::notify(vm);
			GetSummary::getSummary(vcf_filename, output_prefix, trio_info_str);

		} else {

			cout << command_info.str() << endl;
		}
	}

	return 0;
}
