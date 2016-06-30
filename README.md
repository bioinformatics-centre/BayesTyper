# BayesTyper

**bayesTyper**: Method for variant graph genotyping based on exact alignment of k-mers.  
**bayesTyperUtils**: Utility tools for pre- and post-processing Variant Call Format (VCF) files for BayesTyper.

Full documentation of all tools in BayesTyperUtils will follow in July 2016. 

### Compilation

##### Requirements
The following software and libraries needs to be installed in order to compile BayesTyper:
* g++ (C++11 support required)
* [CMake](https://cmake.org/) (version 2.8.0 or higher)
* [Boost](http://www.boost.org) (works with version 1.56.0)

##### BayesTyper
To compile BayesTyper:
1. `cd bayesTyper`
2. `mkdir build`
3. `cd build`
4. `cmake ../src/`
5. `make`

This will create the binary: `bayesTyper/build/bayesTyper`

##### BayesTyperUtils
To compile BayesTyperUtils do the following steps:
1. `cd bayesTyperUtils`
2. `cmake .`
3. `make`

This will create the binary: `bayesTyperUtils/bin/BayesTyperUtils`

### Usage
To run BayesTyper use the following command-line:

`bayesTyper -v <variants> -s <samples> -g <genome> -d <decoy> -p <threads>` 


##### Set of variants (`<variants>`)
Set of variants (VCF format) to be genotyped by BayesTyper (required). Only the 5 first columns in the VCF file are required.  


##### Sample information (`<samples>`)
Tab-delimited text file with sample information, where each line corresponds to a sample (required). The file should contain three columns: 

1. Sample name (only used in output header).
2. Sample gender (either *F* for female or *M* for male).
3. Path to sample k-mer counts ([KMC2](http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=kmc&subpage=about) output). The suffixes *.kme_pre* and *.kmc_suf* should not be included.

Example:

*Sample1&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;F&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\<path_to_sample1_kmers>*  
*Sample2&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;M&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\<path_to_sample2_kmers>*

The k-mer counts should be generated using [KMC2](http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=kmc&subpage=about) version 2.3.0 (see example pipeline).   


##### Reference genome (`<genome>`)
Reference genome (fasta format) containing autosomal and allosomal chromosomes with known ploidy (required).  

##### Decoy sequences (`<decoy>`)
Decoy sequences (fasta format) with unknown ploidy such as the mitochondrial sequence (recommended).   

##### Computational requirements
BayesTyper is multi-threaded and the number of CPU threads can be specified using `-p`. The memory footprint is currently quite high (~15M variants and 10x coverage: **~195GB** for 1 sample and **~365GB** for 10 samples). We expect a marked reduction in memory footprint in the next release that is due in July 2016.

### Example pipeline
The following is an example of a pipeline for genotyping a single sample on variants from two sources (chr22 only). The data can be downloaded [here](http://people.binf.ku.dk/~lassemaretty/bt_example_chr22.tar.gz).

Notice that the pipeline requires **16GB** of memory and 24 threads. The latter can be changed in the [KMC2](http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=kmc&subpage=about) and BayesTyper command-lines.

##### Merging
Merge the two input variant sets (1000 genomes CEU variants and Platypus predictions) into a single call-set:

`./bayesTyperUtils/bin/BayesTyperUtils combine -v CEU:variants_ceu_chr22.vcf,PP:variants_pp_chr22.vcf -o variants_merged_chr22`

##### K-mer counting
Download and install [KMC2](http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=kmc&subpage=about) version 2.3.0. Count all k-mers (k=55) in the sequencing reads:

`kmc -k55 -m16 -sm -t24 -ci1 @reads_chr22.txt kmers_chr22 .`

*reads.txt* contains a list of the fastq files.

##### Genotyping
Genotype the merged variant set using the counted k-mers:

`./bayesTyper/build/bayesTyper -v variants_merged_chr22.vcf -s sample.txt -g chr22.fa -d decoy.fa -p 24` 

*sample.txt* contains the sample name, gender and a path to the k-mer counts.

##### Filtering
Filter genotyped alleles with insufficient k-mer information:

`./bayesTyperUtils/bin/BayesTyperUtils filter -v bayesTyper_genotypes.vcf -o bayesTyper_genotypes_filtered --min-nok 3 --min-gpp 0.9 --filter-dependencies --keep-filtered-variants` 








