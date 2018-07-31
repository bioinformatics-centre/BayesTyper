# BayesTyper #
BayesTyper performs genotyping of all types of variation (including SNPs, indels and complex structural variants) based on an input set of variants and read k-mer counts. Internally, BayesTyper uses exact alignment of k-mers to a graph representation of the input variants and reference sequence in combination with a probabilistic model of k-mer counts to do genotyping.

## News ##
* 31 July 2018: BayesTyper data bundles have been updated to fix reference naming issue when used with the snakemake workflow and to fix a fasta parsing error in `bayesTyperTools convertAllele`
* 19 June 2018: BayesTyper *canon* and *decoy* references in the data bundles have been updated
* 18 June 2018: BayesTyper paper published as *Technical Report* in *Nature Genetics* ([link](https://www.nature.com/articles/s41588-018-0145-5)). 
    * Please note that the results presented in the paper were obtained using BayesTyper [v1.2.0](https://github.com/bioinformatics-centre/BayesTyper/releases/tag/v1.2) 
    * The new BayesTyper [v >= 1.3](https://github.com/bioinformatics-centre/BayesTyper/releases/latest) requires **much less memory** and **runs faster** than stated in *Supplementary Table 2* in the manuscript with no impact on sensitivity and genotyping acceracy as compared with v1.2. Please refer to the *Computational requirements* section below and the [release notes](https://github.com/bioinformatics-centre/BayesTyper/releases/tag/v1.3) for version 1.3 for further details. 

* 18 June 2018: Patch ([v1.3.1](https://github.com/bioinformatics-centre/BayesTyper/releases/tag/v1.3.1)) 
    
* 15 June 2018: New major release out ([v1.3](https://github.com/bioinformatics-centre/BayesTyper/releases/tag/v1.3)) featuring:
    * **New interface:** bayesTyper has been refactored into `BayesTyper cluster` and `BayesTyper genotype`. The *cluster* part partitions the variants into *units* that can then be *genotyped*.
    * **Much reduced memory:** Only graphs and k-mers for a single unit need to reside in memory at the same time; the rest remains on disk. A bloom filter stores information about k-mers shared across units.  This construct ensures that *memory usage is (almost) independent of the number of candidate variants*.
    * **Cluster support:** Each unit can be genotyped independently and hence distributed across nodes on a cluster followed by simple concatenation of the unit vcf files (e.g. using `bcftools concat`).
    * **Simultaneous genotyping and filtering:** No need to run `bayesTyperTools filter`. Hard filters are now applied up front by `bayesTyper genotype`. Genotypes can still be refiltered using `bayesTyperTools filter` after genotyping if necessary.
    * **Parallel read bloom generation:** `bayesTyperTools makeBloom` can now use multiple threads (and  scales very well with the number of threads).
    * **snakemake workflow:** We have added an example *snakemake* workflow to the repo - this can orchestrate the entire pipeline straight from BAM(s) over variant candidates to final genotypes.  
    
* 26 February 2018: Manuscript release ([v1.2](https://github.com/bioinformatics-centre/BayesTyper/releases/tag/v1.2))

## Why should I use BayesTyper? ##
### Short explanation ###
Because it allows you to obtain accurate genotypes spanning from SNVs/short indels to complex structural variants and hence provides a more complete picture of the genome as compared with standard methods - without sacrificing accuracy.

### Detailed explanation ###
Standard methods for genotyping (e.g. *GATK-HaplotypeCaller*, *Platypus* and *Freebayes*) start from an alignment of reads (e.g. by *BWA-MEM*) and then
1. Call candidate variants.
2. Perform *local* realignment of reads anchored to a particular variable region to candidate variant haplotypes.
3. Estimate genotypes (and thus final variant calls).

This approach can result in a bias towards the reference sequence since reads informative for a particular variant may have been left either unaligned (because of too large an edit distance to the reference) or have aligned better elsewhere in the reference.

The variant graph approach used by BayesTyper ensures that the resulting calls are not biased towards the reference sequence by effectively realigning *all* reads (or more specifically their k-mers) when genotyping candidate variants. In our recent [paper](https://www.nature.com/articles/s41588-018-0145-5), we show how this approach significantly improves both sensitivity and genotyping accuracy for most variant types - especially non-SNVs (please see citation below).

## Installation ##
1. Download the latest static Linux x86_64 build (k=55) found under [releases](https://github.com/bioinformatics-centre/BayesTyper/releases/latest).
    * To build from source, please refer to the [build wiki](https://github.com/bioinformatics-centre/BayesTyper/wiki/Building-BayesTyper-from-source) for detailed build instructions. Note that the k-mer size is determined compile time. Hence, if you want k ≠ 55 you need to compile it yourself (or post a [feature request](https://github.com/bioinformatics-centre/BayesTyper/issues) and we will see what we can do).

2. Download the BayesTyper human data bundle ([GRCh37](http://people.binf.ku.dk/~lassemaretty/bayesTyper/bayestyper_GRCh37_bundle.tar.gz) and [GRCh38](http://people.binf.ku.dk/~lassemaretty/bayesTyper/bayestyper_GRCh38_bundle.tar.gz)) containing reference sequences preprocessed for BayesTyper (i.e. canonical and decoy chromosomes) together with a reference matched [variation prior](https://github.com/bioinformatics-centre/BayesTyper/wiki/Variation-prior) database.

## Usage ##
The BayesTyper genotyping process occurs in two stages:
1. Generation of variant candidates (using other tools)
2. Genotyping based on variant candidates (using *BayesTyper*)

As indicated, BayesTyper **does not** find candidate variants on its own. Instead, users can combine the variant discovery strategies suitable for their study as it will depend on the study design (e.g. coverage, number of samples etc.) as well as the available resources.

Below we outline an example strategy, where candiates are obtained using
* *GATK-HaplotypeCaller* (without all the pre-processing of alignments and post-processing of variants)
* *Platypus*
* *Manta*
* The [variation prior](https://github.com/bioinformatics-centre/BayesTyper/wiki/Variation-prior) (part of the BayesTyper data bundle, see [Installation](https://github.com/bioinformatics-centre/BayesTyper#installation))

The complete workflow (i.e. BAM(s) to genotypes) outlined below is further provided as a [snakemake workflow](https://github.com/bioinformatics-centre/BayesTyper/tree/master/workflows) for easy deployment of BayesTyper. Please refer to the [snakemake wiki](https://github.com/bioinformatics-centre/BayesTyper/wiki/Running-BayesTyper-using-snakemake) for detailed instructions on how to set up and execute the workflow on your data.

**Important:** This workflow should work well for most cases. If you prefer to use another approach, please note that the candidate variant set *must contain* at least 1 million *candidate* SNVs that are needed for accurate estimation of parameters.

**Important:** Please note that it is currently only possible to genotype 30 samples at the time using *BayesTyper*. To run more samples, please execute *BayesTyper* in batches as described in the [batching wiki](https://github.com/bioinformatics-centre/BayesTyper/wiki/Executing-BayesTyper-on-sample-batches). Batching is currently not supported by the *snakemake* workflow - please let us know if you require this feature by filing a [feature request](https://github.com/bioinformatics-centre/BayesTyper/issues).

**Important:** Bayestyper supports uncompressed and gzip compressed vcf files. Please note that bgzip compression is currently **not** supported.

**Important:** If you intend to run on **non-human** data, please refer to the [Running on non-human data wiki](https://github.com/bioinformatics-centre/BayesTyper/wiki/Running-BayesTyper-on-non-human-data) for more information. 

### 1. Generation of variant candidates ###
Starting from a set of indexed, aligned reads (obtained e.g. using *BWA-MEM*):

1. For each sample, run [*HaplotypeCaller*](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php) to get standard mapping-based candidates
    * Running Base Quality Recalibration or doing joint genotyping is **not** required
    * Faster alternative: [*Freebayes*](https://github.com/ekg/freebayes)
    
2. For each sample, run [*Platypus*](http://www.well.ox.ac.uk/platypus) to identify small and medium sized variants

3. For each sample, run [*Manta*](https://github.com/Illumina/manta) to identify candidates by *de novo* local assembly (important for detecting larger deletions and insertions). Convert allele IDs (e.g. \<DEL\>) in the Manta output to sequences using `bayesTyperTools convertAllele`

4. For each caller, left-align and normalize variants using `bcftools norm` ([*bcftools*](https://github.com/samtools/bcftools))

5. Combine variants across *all* samples, callers and the [variation prior](https://github.com/bioinformatics-centre/BayesTyper/wiki/Variation-prior) using `bayesTyperTools combine -v GATK:<gatk_sample1>.vcf,GATK:<gatk_sample2>.vcf,PLATYPUS:<platypus_sample1>.vcf,PLATYPUS:<platypus_sample2>.vcf,MANTA:<manta_sample1>.vcf,...,prior:<prior>.vcf -o <candiate_variants_prefix> -z`
   * **Note:** The source tag before the colon (e.g. GATK) only serves to identify the origin of the calls in the final callset - it can take any value.
   * **Important:** `bayesTyperTools combine` requires the vcf header to contain contig entries (e.g.`##contig=<ID=8,length=146364022>`) for all reference sequences containing variants in the vcf; the contigs further need to appear in the same order in the header and for the variant entries.

### 2. Genotyping based on variant candidates ###

1. Count k-mers
   1. Run [KMC3](https://github.com/refresh-bio/KMC) on each sample (set k=55 using `-k55` and include singleton k-mers using `-ci1`)
      * **Note:** Default is fq(.gz) input - bam input is enabled using `-fbam`.
      * **Important:** If the reads are in bam format, make sure the file also contains the unmapped reads.   
   2. Create a read k-mer bloom filter for each sample from the KMC3 output using `bayesTyperTools makeBloom -k <kmc_output_prefix> -p <num_threads>`
      * **Important:** The resulting bloom filter (*<sample_id>.bloomMeta* and *<sample_id>.bloomData*) and the KMC3 output (*<sample_id>.kmc_pre* and *<sample_id>.kmc_suf*) must reside in the same directory and have the same prefix

2. Identify variant clusters: `bayesTyper cluster -v <candiate_variants_prefix>.vcf.gz -s <samples>.tsv -g <ref_build>_canon.fa -d <ref_build>_decoy.fa -p <num_threads>`
      * **Note:** This partitions the candidate variants into units written to separate directories (*bayestyper_unit_1*, *bayestyper_unit_2*, ...) - each containing between 5M and 10M variants. These can be genotyped independently e.g. on a cluster (supported by the [snakemake workflow](https://github.com/bioinformatics-centre/BayesTyper/tree/master/workflows)) followed by a simple concatenation operation using `bcftools concat` (see below).
      * **Important:** The `<samples>.tsv` file should contain one sample per row with columns *<sample_id>, \<sex> and <kmc_output_prefix>* and no header ([example](http://people.binf.ku.dk/~lassemaretty/bayesTyper/bt_samples_example.tsv)).
      * **Important:** Data that is common to all units and necessary for genotyping is written to the *bayestyper_cluster_data* directory.
      * **Note:** Human reference files (canon/decoy) are provided in the BayesTyper data bundle (see [Installation](https://github.com/bioinformatics-centre/BayesTyper#installation)).

3. Genotype variant clusters: `bayesTyper genotype -v bayestyper_unit_<unit_id>/variant_clusters.bin -c bayestyper_cluster_data -s <samples>.tsv -g <ref_build>_canon.fa -d <ref_build>_decoy.fa -o bayestyper_unit_<unit_id>/bayestyper -z -p <num_threads>`
      * **Note:** The genotype command also applies the default BayesTyper hard-filters by setting the variant FILTER status and the sample specific allele filter (SAF) format attribute. Please refer to the [filter wiki](https://github.com/bioinformatics-centre/BayesTyper/wiki/Filtering) for information about the filters used, how to changes the defaults up front and how to refilter the data after running the genotyping step.
      * **Important:** The filtering procedure only filters the genotypes ("./.") and hence, downstream tools **should prefilter on the variant quality and FILTER==\"PASS\"** (e.g. using [*bcftools*](https://github.com/samtools/bcftools)) to obtain the filtered calls.

4. Concatenate units using [*bcftools*](https://github.com/samtools/bcftools): `bcftools concat -O z -o <output_prefix>.vcf.gz bayestyper_unit_1/bayestyper.vcf.gz bayestyper_unit_2/bayestyper.vcf.gz ...`
    * **Important:** Unit arguments to bcftools should be in ascending order (unit_1 unit_2 ...) for the output to be properly sorted

## Computational requirements ##

|Number of samples|Coverage|Number of variant alleles|Max allele length (nts)|Number of threads|Wall time (h, single node)|Wall time (h, multiple nodes)\*|Max memory (GB)|
|-|---|-----|-------|--|-|-|--|
|3|13x|21.4M|500,000|28|5-6|2-3|41|
|3|13x|64.4M|500,000|28|17-18|4-5|42|
|13|50x|11.7M|10,000|28|31-32|16-17|66|
|13|50x|61.1M|10,000|28|92-93|15-16|62|

\* `bayesTyper genotype` can be distributed across nodes on a cluster - between 2 and 11 nodes were used in this benchmark.

The time estimates are for running `bayesTyper cluster` and `bayesTyper genotype` only. Expect <1h combined run-time per sample for counting k-mers by *KMC3* and bloom filter creation by *bayesTyperTools*. All runs were done on a 64-bit Intel Xeon 2.00 GHz machine with 128 GB of memory.

## Citing BayesTyper ##
Sibbesen JA*, Maretty L*, The Danish Pan-Genome Consortium & Krogh A: Accurate genotyping accross variant classes and lengths using variant graphs. *Nature Genetics*, 2018. [link](https://www.nature.com/articles/s41588-018-0145-5). *Equal contributors.

## Studies that have used BayesTyper ##
* Sequencing and de novo assembly of 150 genomes from Denmark as a population reference. *Nature*, 2017 ([link](https://www.nature.com/articles/nature23264))
* A high-coverage Neandertal genome from Vindija Cave in Croatia. *Science*, 2017 ([link](http://science.sciencemag.org/content/early/2017/10/04/science.aao1887))
* Analysis of 62 hybrid assembled human Y chromosomes exposes rapid structural changes and high rates of gene conversion. *PLOS Genetics*, 2017 ([link](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006834))
* Assembly and analysis of 100 full MHC haplotypes from the Danish population, *Genome Research*, 2017 ([link](https://genome.cshlp.org/content/early/2017/08/03/gr.218891.116))

Please let us know if you use BayesTyper in your publication - then we will put it on the list.

## Contact ##
Please post an [issue](https://github.com/bioinformatics-centre/BayesTyper/issues) if you have questions regarding how to run BayesTyper, if you want to report bugs or request new features. You can also reach us at jasi at binf dot ku dot dk or lasse dot maretty at clin dot au dot dk.

## Third-party software acknowledgements ##
We thank the developers of the third-party libraries used by BayesTyper:
* [Edlib](https://github.com/Martinsos/edlib)
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
* [KMC](https://github.com/refresh-bio/KMC)
* [ntHash](https://github.com/bcgsc/ntHash)
