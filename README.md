# BayesTyper #
BayesTyper performs genotyping of all types of variation (including structural and complex variation) based on an input set of variants and read k-mer counts. Internally, BayesTyper uses exact alignment of k-mers to a graph representation of the input variation and reference sequence in combination with a probabilistic model of k-mer counts to do genotyping. The variant representation ensures that the resulting calls are not biased towards the reference sequence. 

The BayesTyper was used to integrate mapping- and assembly-based calls in the [GenomeDenmark project](http://www.nature.com/nature/journal/vaop/ncurrent/full/nature23264.html). A manuscript describing the method is currently in revision.

## Use cases ##

### Variant integration ###
Sensitive calling of structural variation typically requires running multiple callers to ensure sensitivity yet this leads to the problem of integrating calls across call-sets. The BayesTyper can be used to produce a fully integrated call-set including SNVs, indels and complex variation from input variant *candidates* produced by a panel of methods; the panel must include standard SNV and indel calls e.g. from GATK, Freebayes or Platypus.

### Prior based genotyping ###
A signficant amount of both simple and complex variation is already known from large population-scale studies. As some of these variants may be missed in a study - even when running multiple methods - due to alignment bias, we provide a database containing common SNPs/indels together with complex variants that can be combined with *in-sample* calls (i.e. calls based only on the study data) to improve sensitivity.
This approach can for instance be used to quickly augment a set of standard SNV and indel calls (e.g. from GATK) with structural variation by running BayesTyper on the SNV/indel calls combined with our variation database. For higher sensitivity, *in-sample* complex variation calls can be combined with the database to produce the final intergrated call-set.

## Installation ##
The BayesTyper package contains *BayesTyper*, which does the genotyping, and *BayesTyperTools*, which is used to pre- and post-process VCF files for BayesTyper.

#### Prerequisites ####
* gcc (c++11 support required. Tested with gcc 4.8 and 4.9 work)
* CMake (version 2.8.0 or higher)
* Boost (tested with version 1.55.0 and 1.56.0)

#### Building BayesTyper ####
BayesTyper currently needs to be build from source; a pre-compiled version will be released at a later time. 
1. `git clone https://github.com/bioinformatics-centre/BayesTyper.git`
2. `cd BayesTyper`
2. `mkdir build && cd build`
5. `cmake ..`
6. `make`

## Basic usage ##
1. Count k-mers for each sample using [kmc3](http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=kmc&subpage=download)
   * e.g. `kmc -k55 sample_1.fq sample_1`
   * For low coverage data (<20X), include singleton k-mers by using `-ci1` instead of `-ci2` 
2. Prepare a tsv file with sample information. One sample per row with columns "sample_id", "sex" and "path_to_kmc3_output" (e.g. sample_1 &nbsp;&nbsp;&nbsp; F &nbsp;&nbsp;&nbsp; kmc3_counts/sample_1)
3. Prepare the variant input
   * `bayesTyperTools combine -o bayesTyper_input -v gatk1:sample_1_gatk.vcf,gatk2:sample_2_gatk.vcf,varDB:SNP_dbSNP150common_SV_1000g_dbSNP150all_GDK_GoNL_GTEx_GRCh38.vcf`
4. Run BayesTyper
   * `bayesTyper -o integrated_calls -s samples.tsv -v bayesTyper_input.vcf -g hg38.fa -p <threads> > bayesTyper_log.txt`
5. Filter output
   1. Get coverage stats for filters: `grep "Estimated" bayesTyper_log.txt | cut -f10,18,21 -d ' ' | tr ' ' '\t' > kmer_coverage_estimates.txt`
   2. Run filtering: `bayesTyperTools filter -o integrated_calls_filtered -v integrated_calls.vcf -g hg38.fa --kmer-coverage-filename kmer_coverage_estimates.txt`

## Variant database ##
The variant databases were constructed by combining
* [dbSNP](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC29783/)
   * Rare SNVs were excluded
* [1000 Genomes Project](https://www.nature.com/nature/journal/v526/n7571/full/nature15394.html)
* [Genome of the Netherlands Project](https://www.nature.com/articles/ncomms12989)
   * GRCh38 version has been lifted from GRCh37 using [crossmap](https://www.ncbi.nlm.nih.gov/pubmed/24351709) 
* [Genotype-Tissue Expression (GTEx) Project](http://www.nature.com/ng/journal/v49/n5/full/ng.3834.html)
   * GRCh38 version has been lifted from GRCh37 using [crossmap](https://www.ncbi.nlm.nih.gov/pubmed/24351709) 
* [GenomeDenmark](http://www.nature.com/nature/journal/vaop/ncurrent/full/nature23264.html) project.
   * GRCh37 version has been lifted from GRCh38 using [crossmap](https://www.ncbi.nlm.nih.gov/pubmed/24351709) 

* [BayesTyper_varDB_GRCh37](http://people.binf.ku.dk/~lassemaretty/bayesTyper/SNP_dbSNP150common_SV_1000g_dbSNP150all_GDK_GoNL_GTEx_GRCh37.vcf)
* [BayesTyper_varDB_GRCh38](http://people.binf.ku.dk/~lassemaretty/bayesTyper/SNP_dbSNP150common_SV_1000g_dbSNP150all_GDK_GoNL_GTEx_GRCh38.vcf)

## Memory requirements ## 
* 50,000,0000 variants, 30X coverage, 10 samples, no singleton k-mers -> 340 Gb
* 50,000,0000 variants, 10X coverage, 10 samples -> 480 Gb
