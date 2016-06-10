
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

#include "vcf++/Utils.hpp"
#include "vcf++/CompareOperators.hpp"

#include "convertAlleleIdToSequence.hpp"
#include "clean.hpp"
#include "remove.hpp"
#include "combine.hpp"
#include "split.hpp"
#include "merge.hpp"
#include "annotate.hpp"
#include "calcFilterStats.hpp"
#include "filter.hpp"
#include "markDeNovo.hpp"
#include "writeIndels.hpp"
#include "addRepeatAnnotation.hpp"
#include "generateRelease.hpp"
#include "getSummary.hpp"
#include "calcTrioConcordance.hpp"
#include "selectValidationVariants.hpp"

namespace po = boost::program_options;

static const string parents_trio_regular_expression_default = "^[0-9]+((-01)|(-02)){1}$";

using namespace std;

int main (int argc, char * const argv[]) {

	std::cout << "\n[" << Utils::getLocalTime() << "] " << "You are using BayesTyper Utilities (" << BTU_VERSION << ")\n" << std::endl;

	stringstream command_info;

	command_info << "Usage BayesTyperUtils <command> [options]" << endl;
	command_info << "Commands:\n" << endl;
	command_info << "convertAlleleIdToSequence\tconvert allele ID's to sequence" << endl;
	command_info << "clean\t\t\t\tremove samples, attributes, missing and redundant alleles" << endl;
	command_info << "remove\t\t\t\tremove region(s)" << endl;
	command_info << "combine\t\t\t\tcombine call-sets (vertical)" << endl;
	command_info << "split\t\t\t\tsplit call-set(s)" << endl;
	command_info << "merge\t\t\t\tmerge call-sets (BayesTyper only)" << endl;
	command_info << "annotate\t\t\tannotate alleles" << endl;
	// command_info << "calcFilterStats\t\t\tcalculate filter statistics (BayesTyper only)" << endl;
	command_info << "filter\t\t\t\tfilter samples (BayesTyper only)" << endl;
	// command_info << "markDeNovo\t\t\tmark de novo events (BayesTyper only)" << endl;
	command_info << "writeIndels\t\t\twrite indel sequence fasta file" << endl;
	command_info << "addRepeatAnnotation\t\tadd RepeatMasker annotation" << endl;
	command_info << "generateRelease\t\t\tgenerate cleaned call-sets for release (BayesTyper only)" << endl;
	command_info << "getSummary\t\t\tget call-set summary (BayesTyper only)" << endl;
	command_info << "calcTrioConcordance\t\tcalculate trio concordance" << endl;
	command_info << "selectValidationVariants\tselect variants for validation (BayesTyper only)" << endl;

	if (argc == 1) {

		cout << command_info.str() << endl;

	} else {

		if (strcmp(argv[1],"convertAlleleIdToSequence") == 0) {

			string vcf_filename;
			string genome_filename;
			string output_prefix;

			string mei_filename = "";
			uint max_allele_length;
			bool keep_imprecise_variants;
			bool is_delly_output;

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
				("max-allele-length", po::value<uint>(&max_allele_length)->default_value(3000000), "maximum allele length (both reference and alternative allele)")
				("keep-imprecise-variants", po::value<bool>(&keep_imprecise_variants)->default_value(false)->implicit_value(true), "do not filter imprecise variants")
				("is-delly-output", po::value<bool>(&is_delly_output)->default_value(false)->implicit_value(true), "vcf file is output from the Delly2 tool.")
			;

			po::options_description desc("## BayesTyperUtils convertAlleleIdToSequence ##");
			desc.add(general).add(required_options).add(optional_options);

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);

			if (vm.count("help") || argc == 2) {
			   cout << desc << endl;
			   return 1;
			}

			po::notify(vm);

			ConvertAlleleIdToSequence::convertAlleleIdToSequence(vcf_filename, genome_filename, output_prefix, mei_filename, max_allele_length, keep_imprecise_variants, is_delly_output);

		} else if (strcmp(argv[1],"clean") == 0) {

			string vcf_filename;
			string output_prefix;

			po::options_description general("", 150);
			general.add_options()
			   ("help,h", "produce help message for base options")
			;

			po::options_description required_options("== Required options ==", 160);
			required_options.add_options()

				("vcf-filename,v", po::value<std::string>(&vcf_filename)->required(), "variant file.")
				("output-prefix,o", po::value<std::string>(&output_prefix)->required(), "output prefix (suffix: \".vcf\").")
			;

			po::options_description desc("## BayesTyperUtils clean ##");
			desc.add(general).add(required_options);

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);

			if (vm.count("help") || argc == 2) {
			   cout << desc << endl;
			   return 1;
			}

			po::notify(vm);

			Clean::clean(vcf_filename, output_prefix);

		} else if (strcmp(argv[1],"remove") == 0) {

			string vcf_filename;
			string region_string;
			string output_prefix;

			po::options_description general("", 150);
			general.add_options()
			   ("help,h", "produce help message for base options")
			;

			po::options_description required_options("== Required options ==", 160);
			required_options.add_options()

				("vcf-filename,v", po::value<std::string>(&vcf_filename)->required(), "variant file.")
				("regions,r", po::value<std::string>(&region_string)->required(), "regions (<chr1>,<start1>,<end1>:<chr2>:...); for whole chromosomes the start and end positions can be omitted.")
				("output-prefix,o", po::value<std::string>(&output_prefix)->required(), "output prefix (suffix: \".vcf\").")
			;

			po::options_description desc("## BayesTyperUtils remove ##");
			desc.add(general).add(required_options);

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);

			if (vm.count("help") || argc == 2) {
			   cout << desc << endl;
			   return 1;
			}

			po::notify(vm);

			Remove::remove(vcf_filename, region_string, output_prefix);

		} else if (strcmp(argv[1],"combine") == 0) {

			string vcf_filenames_str;
			string output_prefix;
			bool do_not_remove_asmvar_snps;

			po::options_description general("", 150);
			general.add_options()
			   ("help,h", "produce help message for base options")
			;

			po::options_description required_options("== Required options ==", 160);
			required_options.add_options()

				("vcf-filenames,v", po::value<std::string>(&vcf_filenames_str)->required(), "call-set names and file (<call-set>:<filename>,<call-set>:<filename>,...). For an annotation use \"Ann\" followed by an optional \"_<name>\" to include identifiers. For AsmVar use \"AsmVar\" followed by an optional \"_<name>\" to include scaffold information (needs genotypes).")
				("output-prefix,o", po::value<std::string>(&output_prefix)->required(), "output prefix (suffix: \".vcf\").")
			;

			po::options_description optional_options("== Optional options ==", 160);
			optional_options.add_options()

				("do-not-remove-asmvar-snps", po::value<bool>(&do_not_remove_asmvar_snps)->default_value(false)->implicit_value(true), "do not remove SNPs from AsmVar call-sets.")
			;

			po::options_description desc("## BayesTyperUtils combine ##");
			desc.add(general).add(required_options).add(optional_options);

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);

			if (vm.count("help") || argc == 2) {
			   cout << desc << endl;
			   return 1;
			}

			po::notify(vm);

			Combine::combine(Utils::splitString(vcf_filenames_str, ','), output_prefix, !do_not_remove_asmvar_snps);

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

				("vcf-filenames,v", po::value<std::string>(&vcf_filenames_str)->required(), "variant files (<call-set>:<filename>,<call-set>:<filename>,...). The number of threads will be equal to the number files.")
			;

			po::options_description optional_options("== Optional options ==", 160);
			optional_options.add_options()

				("min-batch-size", po::value<uint>(&min_batch_size)->default_value(50000), "minimum number of variants in a batch.")
				("use-vcg-split", po::value<bool>(&use_vcg_split)->default_value(false)->implicit_value(true), "require split between variant cluster groups.")
				("use-chr-split", po::value<bool>(&use_chr_split)->default_value(false)->implicit_value(true), "require split between chromosomes (overwrites <min-batch-size>).")
			;

			po::options_description desc("## BayesTyperUtils split ##");
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
			string parents_trio_regular_expression;

			po::options_description general("", 150);
			general.add_options()
			   ("help,h", "produce help message for base options")
			;

			po::options_description required_options("== Required options ==", 160);
			required_options.add_options()

				("vcf-filenames,v", po::value<std::string>(&vcf_filenames_str)->required(), "variant files.")
				("output-prefix,o", po::value<std::string>(&output_prefix)->required(), "output prefix (suffix: \".vcf\").")
			;

			po::options_description desc("## BayesTyperUtils merge ##");
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
			float min_allele_overlap;
			bool no_sequence_match;

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

				("min-allele-overlap", po::value<float>(&min_allele_overlap)->default_value(0.5, "0.5"), "minimum overlap between annotation allele and predicted allele.")
				("no-sequence-match", po::value<bool>(&no_sequence_match)->default_value(false)->implicit_value(true), "do not require sequence match.")
			;

			po::options_description desc("## BayesTyperUtils annotate ##");
			desc.add(general).add(required_options).add(optional_options);

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);

			if (vm.count("help") || argc == 2) {
			   cout << desc << endl;
			   return 1;
			}

			po::notify(vm);

			Annotate::annotate(vcf_filename, annotation_filename, output_prefix, min_allele_overlap, !no_sequence_match);

		} else if (strcmp(argv[1],"calcFilterStats") == 0) {

			cerr << "\nERROR: Incompatible with the latest version of BayesTyper.\n" << endl;
			return 1;	
			
			// string vcf_filename;
			// string output_prefix;

			// string nok_thresholds;
			// string nuk_thresholds;

			// string parents_trio_regular_expression;
			// float min_called_probability;

			// po::options_description general("", 150);
			// general.add_options()
			//    ("help,h", "produce help message for base options")
			// ;

			// po::options_description required_options("== Required options ==", 160);
			// required_options.add_options()

			// 	("vcf-filename,v", po::value<std::string>(&vcf_filename)->required(), "variant file.")
			// 	("output-prefix,o", po::value<std::string>(&output_prefix)->required(), "output prefix (suffix: \"_stats.txt\").")
			// ;

			// po::options_description optional_options("== Optional options ==", 160);
			// optional_options.add_options()

			// 	("nok-thresholds", po::value<string>(&nok_thresholds)->default_value("0"), "filter genotypes were at least one allele has less than <threshold> observed (count above zero) kmers (comma separated).")
			// 	("nuk-thresholds", po::value<string>(&nuk_thresholds)->default_value("0"), "filter genotypes were at least one allele has less than <threshold> unique (not intercluster or multicluster) kmers (comma separated).")
			// 	("parents-trio-reg-expression", po::value<string>(&parents_trio_regular_expression)->default_value(parents_trio_regular_expression_default), "regular expression for matching parents in trio samples (used when calculating population statistics).")
			// 	("min-call-probability", po::value<float>(&min_called_probability)->default_value(0.9, "0.9"), "minimum probability needed for an allele or genotype to be considered called.")
			// ;

			// po::options_description desc("## BayesTyperUtils calcFilterStats ##");
			// desc.add(general).add(required_options).add(optional_options);

			// po::variables_map vm;
			// po::store(po::parse_command_line(argc, argv, desc), vm);

			// if (vm.count("help") || argc == 2) {
			//    cout << desc << endl;
			//    return 1;
			// }

			// po::notify(vm);

			// CalcFilterStats::calcFilterStats(vcf_filename, output_prefix, parents_trio_regular_expression, min_called_probability, Utils::splitString(nk_thresholds, ','), Utils::splitString(nok_thresholds, ','), Utils::splitString(nuk_thresholds, ','));

		} else if (strcmp(argv[1],"filter") == 0) {

			string vcf_filename_str;
			string output_prefix;

			float min_nok;
			float min_nuk;
			float min_gpp;

			bool filter_dependencies;
			bool keep_filtered_variants;

			po::options_description general("", 150);
			general.add_options()
			   ("help,h", "produce help message for base options")
			;

			po::options_description required_options("== Required options ==", 160);
			required_options.add_options()

				("vcf-filename,v", po::value<std::string>(&vcf_filename_str)->required(), "variant file.")
				("output-prefix,o", po::value<std::string>(&output_prefix)->required(), "output prefix (suffix: \".vcf\").")
			;

			po::options_description optional_options("== Optional options ==", 160);
			optional_options.add_options()

				("min-nok", po::value<float>(&min_nok), "filter sampled alleles with less than <value> observed (count above zero) kmers.")
				("min-nuk", po::value<float>(&min_nuk), "filter sampled alleles with less than <value> unique (not intercluster or multicluster) kmers.")
				("min-gpp", po::value<float>(&min_gpp)->default_value(0.90, "0.90"), "filter genotypes with a maximum genotype posterior probability lower than <value>.")
				("filter-dependencies", po::value<bool>(&filter_dependencies)->default_value(false)->implicit_value(true), "filter genotypes if one of its dependencies are filtered (nested variants).")
				("keep-filtered-variants", po::value<bool>(&keep_filtered_variants)->default_value(false)->implicit_value(true), "keep filtered variants.")
			;

			po::options_description desc("## BayesTyperUtils filter ##");
			desc.add(general).add(required_options).add(optional_options);

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);

			if (vm.count("help") || argc == 2) {
			   cout << desc << endl;
			   return 1;
			}

			po::notify(vm);

			vector<pair<AttributeFilter*, AttributeFilter*> > attribute_filters;

			Attribute::IsGreaterOrEqCmpOp is_greater_or_eq;
			Attribute::IsEqCmpOp is_eq;

			float missing_value = -1;

			if (vm.count("min-nok")) {

				cout << "[" << Utils::getLocalTime() << "] Adding NOK filter: Sampled allele should have at least " << min_nok << " observed (count above zero) kmers." << endl;

				assert(min_nok >= 0);
				attribute_filters.push_back(make_pair(new ValueAttributeFilter("NOK", Attribute::Value(missing_value), &is_eq), new ValueAttributeFilter("NOK", Attribute::Value(min_nok), &is_greater_or_eq)));
			}

			if (vm.count("min-nuk")) {

				cout << "[" << Utils::getLocalTime() << "] Adding NUK filter: Sampled allele should have at least " << min_nuk << " observed (count above zero) kmers." << endl;

				assert(min_nuk >= 0);
				attribute_filters.push_back(make_pair(new ValueAttributeFilter("NUK", Attribute::Value(missing_value), &is_eq), new ValueAttributeFilter("NUK", Attribute::Value(min_nuk), &is_greater_or_eq)));
			}

			cout << "[" << Utils::getLocalTime() << "] Adding GPP filter: Genotype should have at least a maximum genotype posterior probability of " << min_gpp << endl;

			Filter::filter(vcf_filename_str, output_prefix, attribute_filters, min_gpp, keep_filtered_variants, filter_dependencies);

			for (auto & attribute_filter: attribute_filters) {

				delete attribute_filter.first;
				delete attribute_filter.second;
			}

		} else if (strcmp(argv[1],"markDeNovo") == 0) {

			cerr << "\nERROR: Incompatible with the latest version of BayesTyper.\n" << endl;
			return 1;			

			// string vcf_filename;
			// string output_prefix;

			// string trio_info_str;
			// bool use_genome_dk_trio_syntax;

			// float de_novo_zero_inflation_threshold;
			// float de_novo_allelic_balance_deviation;

			// po::options_description general("", 150);
			// general.add_options()
			//    ("help,h", "produce help message for base options")
			// ;

			// po::options_description required_options("== Required options ==", 160);
			// required_options.add_options()

			// 	("vcf-filename,v", po::value<std::string>(&vcf_filename)->required(), "variant file.")
			// 	("output-prefix,o", po::value<std::string>(&output_prefix)->required(), "output prefix (suffix: \".vcf\" and \"_coverage_stats.txt\").")
			// ;

			// po::options_description optional_options("== Optional options ==", 160);
			// optional_options.add_options()

			// 	("trio-info", po::value<std::string>(&trio_info_str), "trio sample id information (<father>,<mother>,<child>:<father>,<mother>,<child>:...). <trio-info> and <use-genome-dk-trio-syntax> are mutually exclusive.")
			// 	("use-genome-dk-trio-syntax", po::value<bool>(&use_genome_dk_trio_syntax)->default_value(false)->implicit_value(true), "use genomeDK trio sample id syntax. <trio-info> and <use-genome-dk-trio-syntax> are mutual exclusive.")
			// 	("de-novo-zero-inflation-threshold", po::value<float>(&de_novo_zero_inflation_threshold)->default_value(0, "0"), "maximum allowed zero inflation allowed in a trio for de novo events.")
			// 	("de-novo-allelic-balance-deviation", po::value<float>(&de_novo_allelic_balance_deviation)->default_value(0.2, "0.2"), "maximum allowed deviation from perfect allelic balance in the child for de novo events (perfect balance = 0.5). Used to differentiate between somatic and germline mutations.")
			// ;

			// po::options_description desc("## BayesTyperUtils markDeNovo ##");
			// desc.add(general).add(required_options).add(optional_options);

			// po::variables_map vm;
			// po::store(po::parse_command_line(argc, argv, desc), vm);

			// if (vm.count("help") || argc == 2) {
			//    cout << desc << endl;
			//    return 1;
			// }

			// po::notify(vm);

			// assert(vm.count("use-genome-dk-trio-syntax"));

			// if (vm.count("trio-info") and use_genome_dk_trio_syntax) {

			// 	cerr << "\nERROR: Options <trio-info> and <use-genome-dk-trio-syntax> are mutually exclusive.\n" << endl;
			// 	return 1;
			// }

			// if (!(vm.count("trio-info")) and !use_genome_dk_trio_syntax) {

			// 	cerr << "\nERROR: No trio information provided (use either the option <trio-info> or <use-genome-dk-trio-syntax>).\n" << endl;
			// 	return 1;
			// }

			// MarkDeNovo::markDeNovo(vcf_filename, output_prefix, trio_info_str, use_genome_dk_trio_syntax, de_novo_zero_inflation_threshold, de_novo_allelic_balance_deviation);

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

			po::options_description desc("## BayesTyperUtils writeIndels ##");
			desc.add(general).add(required_options).add(optional_options);

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);

			if (vm.count("help") || argc == 2) {
			   cout << desc << endl;
			   return 1;
			}

			po::notify(vm);
			WriteIndels::writeIndels(vcf_filename, output_prefix, min_indel_length);

		} else if (strcmp(argv[1],"addRepeatAnnotation") == 0) {

			string vcf_filename;
			string rm_filename;
			string output_prefix;

			po::options_description general("", 150);
			general.add_options()
			   ("help,h", "produce help message for base options")
			;

			po::options_description required_options("== Required options ==", 160);
			required_options.add_options()

				("vcf-filename,v", po::value<std::string>(&vcf_filename)->required(), "variant file.")
				("rm-file,r", po::value<std::string>(&rm_filename)->required(), "RepeatMasker .out file.")
				("output-prefix,o", po::value<std::string>(&output_prefix)->required(), "output prefix (suffix: \".vcf\")")
			;

			po::options_description desc("## BayesTyperUtils addRepeatAnnotation ##");
			desc.add(general).add(required_options);

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);

			if (vm.count("help") || argc == 2) {
			   cout << desc << endl;
			   return 1;
			}

			po::notify(vm);

			AddRepeatAnnotation::addRepeatAnnotation(vcf_filename, rm_filename, output_prefix);

		} else if (strcmp(argv[1],"generateRelease") == 0) {

			string vcf_filename;
			string output_prefix;
			float min_called_probability;

			po::options_description general("", 150);
			general.add_options()
			   ("help,h", "produce help message for base options")
			;

			po::options_description required_options("== Required options ==", 160);
			required_options.add_options()

				("vcf-filename,v", po::value<std::string>(&vcf_filename)->required(), "variant file.")
				("output-prefix,o", po::value<std::string>(&output_prefix)->required(), "output prefix (suffix: \"_complete.vcf\", \"_called_alleles.vcf\" & \"_de_novo.vcf\").")
			;

			po::options_description optional_options("== Optional options ==", 160);
			optional_options.add_options()

				("min-call-probability", po::value<float>(&min_called_probability)->default_value(0.9, "0.9"), "minimum probability needed for an allele or genotype to be considered called.")
			;

			po::options_description desc("## BayesTyperUtils generateRelease ##");
			desc.add(general).add(required_options).add(optional_options);

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);

			if (vm.count("help") || argc == 2) {
			   cout << desc << endl;
			   return 1;
			}

			po::notify(vm);
			GenerateRelease::generateRelease(vcf_filename, output_prefix, min_called_probability);

		} else if (strcmp(argv[1],"getSummary") == 0) {

			string vcf_filename;
			string output_prefix;
			float min_called_probability;
			string parents_trio_regular_expression;
			string excluded_sample_ids;

			po::options_description general("", 150);
			general.add_options()
			   ("help,h", "produce help message for base options")
			;

			po::options_description required_options("== Required options ==", 160);
			required_options.add_options()

				("vcf-filename,v", po::value<std::string>(&vcf_filename)->required(), "variant file.")
				("output-prefix,o", po::value<std::string>(&output_prefix)->required(), "output prefix (suffix: \"_variant.txt\", \"_allele.txt\", \"_mac.txt\" & \"_de_novo.txt\").")
			;

			po::options_description optional_options("== Optional options ==", 160);
			optional_options.add_options()

				("min-call-probability", po::value<float>(&min_called_probability)->default_value(0.9, "0.9"), "minimum probability needed for an allele or genotype to be considered called.")
				("parents-trio-reg-expression", po::value<string>(&parents_trio_regular_expression)->default_value(parents_trio_regular_expression_default), "regular expression for matching parents in trio samples (used when calculating population statistics).")
				("excluded-sample-ids", po::value<string>(&excluded_sample_ids)->default_value(""), "comma seperated list of samples not to be included in the summary output.")
			;

			po::options_description desc("## BayesTyperUtils getSummary ##");
			desc.add(general).add(required_options).add(optional_options);

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);

			if (vm.count("help") || argc == 2) {
			   cout << desc << endl;
			   return 1;
			}

			po::notify(vm);
			GetSummary::getSummary(vcf_filename, output_prefix, min_called_probability, parents_trio_regular_expression, Utils::splitString(excluded_sample_ids, ','));

		} else if (strcmp(argv[1],"calcTrioConcordance") == 0) {

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
				("output-prefix,o", po::value<std::string>(&output_prefix)->required(), "output prefix (suffix: \"_concordance_trio_stats.txt\").")
			;

			po::options_description optional_options("== Optional options ==", 160);
			optional_options.add_options()

				("trio-info", po::value<std::string>(&trio_info_str), "trio sample id information (<father>,<mother>,<child>:<father>,<mother>,<child>:...). If not specified the following syntax is assumed: father = <name>-01, mother = <name>-02 & child = <name>-(03-99).")
			;

			po::options_description desc("## BayesTyperUtils calcTrioConcordance ##");
			desc.add(general).add(required_options).add(optional_options);

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);

			if (vm.count("help") || argc == 2) {
			   cout << desc << endl;
			   return 1;
			}

			po::notify(vm);

			CalcTrioConcordance::calcTrioConcordance(vcf_filename, output_prefix, trio_info_str);

		} else if (strcmp(argv[1],"selectValidationVariants") == 0) {

			string vcf_filename;
			string val_trio_ids_str;
			string genome_fasta_filename;
			string output_prefix;
			uint min_flank_length;
			uint max_fragment_length;
			float min_gpp;
			float min_acp;

			po::options_description general("", 150);
			general.add_options()
			   ("help,h", "produce help message for base options")
			;

			po::options_description required_options("== Required options ==", 160);
			required_options.add_options()

				("vcf-filename,v", po::value<std::string>(&vcf_filename)->required(), "variant file.")
				("val-trio-ids,s", po::value<std::string>(&val_trio_ids_str)->required(), "comma separated list of ids for trios used for validation. Use \".\" to indicate that all trios should be used.")
				("genome-fasta-filename,g", po::value<std::string>(&genome_fasta_filename)->required(), "genome fasta file")
				("output-prefix,o", po::value<std::string>(&output_prefix)->required(), "output prefix")
			;

			po::options_description optional_options("== Optional options ==", 160);
			optional_options.add_options()

				("min-flank-length", po::value<uint>(&min_flank_length)->default_value(25), "minimum free flank of each side of variant")
				("max-fragment-length", po::value<uint>(&max_fragment_length)->default_value(2000), "maximum length of PCR fragment")
				("min-gpp", po::value<float>(&min_gpp)->default_value(0.9, "0.9"), "mininum genotype posterior probability required.")
				("min-acp", po::value<float>(&min_acp)->default_value(0.9, "0.9"), "mininum allele call posterior probability required.")
			;

			po::options_description desc("## BayesTyperUtils selectValidationVariants ##");
			desc.add(general).add(required_options).add(optional_options);

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);

			if (vm.count("help") || argc == 2) {
			   cout << desc << endl;
			   return 1;
			}

			po::notify(vm);

			vector<string> val_trio_ids;

			if (val_trio_ids_str != ".") {

				val_trio_ids = Utils::splitString(val_trio_ids_str, ',');
				assert(!val_trio_ids.empty());
			}

			SelectValidationVariants::selectValidationVariants(vcf_filename, val_trio_ids, genome_fasta_filename, output_prefix, min_flank_length, max_fragment_length, min_gpp, min_acp);

		} else {

			cout << command_info.str() << endl;
		}
	}

	return 0;
}
