# BayesTyper #
BayesTyper performs genotyping of all types of variation (including structural and complex variation) based on an input set of variants and read k-mer counts. Internally, BayesTyper uses exact alignment of k-mers to a graph representation of the input variation and reference sequence in combination with a probabilistic model of k-mer counts to do genotyping. The variant representation ensures that the resulting calls are not biased towards the reference sequence. 

The BayesTyper was used to integrate mapping- and assembly-based calls in the [GenomeDenmark project](http://www.nature.com/nature/journal/vaop/ncurrent/full/nature23264.html). A manuscript describing the method is currently in revision.

The BayesTyper is being developed by Jonas Andreas Sibbesen, Lasse Maretty and Anders Krogh at the Section for Computational and RNA Biology, Department of Biology, University of Copenhagen.

## Use cases ##

### Variant integration ###
Sensitive calling of structural variation typically requires running multiple callers to ensure sensitivity yet this leads to the problem of integrating calls across call-sets. The BayesTyper can be used to produce a fully integrated call-set including SNVs, indels and complex variation from input variant *candidates* produced by a panel of methods; the panel must include standard SNV and indel calls e.g. from GATK, Freebayes or Platypus.

### Prior based genotyping ###
A signficant amount of both simple and complex variation is already known from large population-scale studies. As some of these variants may be missed in a study - even when running multiple methods - due to alignment bias, we provide a database containing common SNPs/indels together with complex variants that can be combined with *in-sample* calls (i.e. calls based only on the study data) to improve sensitivity.
This approach can for instance be used to quickly augment a set of standard SNV and indel calls (e.g. from GATK) with structural variation by running BayesTyper on the SNV/indel calls combined with our variation database. For higher sensitivity, *in-sample* complex variation calls can be combined with the database to produce the final intergrated call-set.

## Installation ##
The BayesTyper package contains `bayesTyper`, which does the genotyping, and `bayesTyperTools`, which is used to pre- and post-process VCF files for BayesTyper.

#### Prerequisites ####
* gcc (c++11 support required. Tested with gcc 4.8 and 4.9)
* CMake (version 2.8.0 or higher)
* Boost (tested with version 1.55.0 and 1.56.0)

#### Building BayesTyper ####
BayesTyper currently needs to be build from source; a pre-compiled version will be released at a later time. 
1. `git clone https://github.com/bioinformatics-centre/BayesTyper.git`
2. `cd BayesTyper`
2. `mkdir build && cd build`
5. `cmake ..`
6. `make`

The compiled `bayesTyper` and `bayesTyperTools` binaries are now located in the `bin` directory.

## Basic usage ##
1. Count k-mers

   1. Run [KMC3](https://github.com/refresh-bio/KMC) on each sample: `kmc -k55 sample_1.fq sample_1` 
      * This will output k-mer counts to `sample_1.kmc_pre` and `sample_1.kmc_suf`.
      * For low coverage data (<20X), include singleton k-mers by adding `-ci1` to the `kmc3` commandline.
      
2. Prepare variant input

   1. Convert allele IDs (e.g. \<DEL\>) to sequence: `bayesTyperTools convertAlleleId -o sample_1_sv_calls_seq -v sample_1_sv_calls.vcf -g hg38.fa`
      * Currently \<DEL\>, \<DUP\>, \<CN[digit(s)]\>, \<CNV\>, \<INV\>, \<INS:ME:[sequence name]\> are supported. The latter require a fasta file with the mobile element insertion sequences.
      * This step can be skipped if the variant sets does not include any allele IDs (e.g. GATK output).
  
   2. Normalise variants using [Bcftools](https://samtools.github.io/bcftools/): `bcftools norm -o sample_1_gatk_norm.vcf -f hg38.fa sample_1_gatk.vcf`
    
   3. Combine variant sets: `bayesTyperTools combine -o bayesTyper_input -v gatk:sample_1_gatk_norm.vcf,gatk:sample_2_gatk_norm.vcf,gatk:sample_3_gatk_norm.vcf,varDB:SNP_dbSNP150common_SV_1000g_dbSNP150all_GDK_GoNL_GTEx_GRCh38.vcf`
      * The contig fields in the headers need to be identical between variant sets and the variants sorted in the same order as the fields.
      
3. Genotype variants

   1. Prepare sample information: Create tsv file with one sample per row with columns \<sample_id\>, \<sex\> and \<path_to_kmc3_output\> ([example](http://people.binf.ku.dk/~lassemaretty/bayesTyper/bt_samples_example.tsv))
   
   2. Run BayesTyper: `bayesTyper -o integrated_calls -s samples.tsv -v bayesTyper_input.vcf -g hg38.fa -p <threads> > bayesTyper_log.txt`
   
4. Filter output

   1. Get coverage estimates for filters: `grep "Estimated" bayesTyper_log.txt | cut -f10,18,21 -d ' ' | tr ' ' '\t' > kmer_coverage_estimates.txt`
   
   2. Run filtering: `bayesTyperTools filter -o integrated_calls_filtered -v integrated_calls.vcf -g hg38.fa --kmer-coverage-filename kmer_coverage_estimates.txt`

## Variant databases ##
* [BayesTyper_varDB_GRCh37](http://people.binf.ku.dk/~lassemaretty/bayesTyper/SNP_dbSNP150common_SV_1000g_dbSNP150all_GDK_GoNL_GTEx_GRCh37.vcf)
* [BayesTyper_varDB_GRCh38](http://people.binf.ku.dk/~lassemaretty/bayesTyper/SNP_dbSNP150common_SV_1000g_dbSNP150all_GDK_GoNL_GTEx_GRCh38.vcf)

### Variant database sources ###
#### GRCh37 ####
|Source|Version|Filters|Lifted|Reference|
|------|-------|-------|------|---------|
|dbSNP|150|No rare SNVs|No|[link](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC29783/)|
|1000 Genomes Project (1KG)|Phase 3|No SNVs|No|[link](https://www.nature.com/nature/journal/v526/n7571/full/nature15394.html)||
|Genome of the Netherlands Project (GoNL)|Release 6|No SNVs|No|[link](https://www.nature.com/articles/ncomms12989)|
|Genotype-Tissue Expression (GTEx) Project|GTEx Analysis V6|No SNVs|No|[link](http://www.nature.com/ng/journal/v49/n5/full/ng.3834.html)|
|GenomeDenmark (GDK)|v1.0|No SNVs|From GRCh38|[link](http://www.nature.com/nature/journal/vaop/ncurrent/full/nature23264.html)|

#### GRCh38 ####
|Source|Version|Filters|Lifted|Reference|
|------|-------|-------|------|---------|
|dbSNP|150|No rare SNVs|No|[link](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC29783/)|
|1000 Genomes Project (1KG)|Phase 3|No SNVs|No|[link](https://www.nature.com/nature/journal/v526/n7571/full/nature15394.html)||
|Genome of the Netherlands Project (GoNL)|Release 6|No SNVs|From GRCh37|[link](https://www.nature.com/articles/ncomms12989)|
|Genotype-Tissue Expression (GTEx) Project|GTEx Analysis V6|No SNVs|From GRCh37|[link](http://www.nature.com/ng/journal/v49/n5/full/ng.3834.html)|
|GenomeDenmark (GDK)|v1.0|No SNVs|No|[link](http://www.nature.com/nature/journal/vaop/ncurrent/full/nature23264.html)|
   
## Memory requirements ## 
|Variants|Coverage|Samples|Singletons included|Threads|Memory (GB)|Time (wall-time hours)|
|--------|--------|-------|-------------------|-------|-----------|----------------------|
|50M|30X|10|No|24|340|67|
|50M|10X|10|Yes|32|480|58|

## Third-party ## 
Third-party software used by BayesTyper.
* [Edlib](https://github.com/Martinsos/edlib)
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
* [KMC](https://github.com/refresh-bio/KMC)
