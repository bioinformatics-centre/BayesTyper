# BayesTyper #
BayesTyper performs genotyping of all types of variation (including SNPs, indels and complex structural variants) based on an input set of variants and read k-mer counts. Internally, BayesTyper uses exact alignment of k-mers to a graph representation of the input variants and reference sequence in combination with a probabilistic model of k-mer counts to do genotyping.

## Installation ##
1. Download the latest static Linux x86_64 build (k=55) found under [releases](https://github.com/bioinformatics-centre/BayesTyper/releases/latest).
    * To build from source, please refer to the [build wiki](https://github.com/bioinformatics-centre/BayesTyper/wiki/Building-BayesTyper-from-source) for detailed build instructions. Note that the k-mer size is determined compile time. Hence, if you want k ≠ 55 you need to compile it yourself (or post a [feature request](https://github.com/bioinformatics-centre/BayesTyper/issues) and we will see what we can do).

2. Download the BayesTyper data bundle ([GRCh37](http://people.binf.ku.dk/~lassemaretty/bayesTyper/bayestyper_GRCh37_bundle.tar.gz) and [GRCh38](http://people.binf.ku.dk/~lassemaretty/bayesTyper/bayestyper_GRCh38_bundle.tar.gz)) containing reference sequences preprocessed for BayesTyper (i.e. canonical and decoy chromosomes) together with a reference matched variation prior database.

## Usage ##
The BayesTyper genotyping process occurs in two stages:
1. Generation of variant candidates (using other tools)
2. Genotyping based on variant candidates (using *BayesTyper*)

As indicated, BayesTyper **does not** find candidate variants on its own. Instead, users can combine the variant discovery strategies suitable for their study as it will depend on the study design (e.g. coverage, number of samples etc.) as well as the available resources.

Below we outline an example strategy, where candiates are obtained using
* *GATK-HaplotypeCaller* (without all the pre-processing of alignments and post-processing of variants)
* *Platypus*
* *Manta*
* The *variation prior* (part of the BayesTyper data bundle, see [Installation](https://github.com/bioinformatics-centre/BayesTyper#installation))

The complete workflow (i.e. BAM(s) to genotypes) outlined below is further provided as a [snakemake workflow](https://github.com/bioinformatics-centre/BayesTyper/tree/master/workflows) for easy deployment of BayesTyper. Please refer to the [snakemake wiki](https://github.com/bioinformatics-centre/BayesTyper/wiki/Running-BayesTyper-using-snakemake) for detailed instructions on how to set up and execute the workflow on your data.

**Important:** This workflow should work well for most cases. If you prefer to use another approach, please note that the candidate variant set *must contain* at least 1 million SNVs that are needed for accurate estimation of parameters.

**Important:** Please note that it is currently only possible to genotype 30 samples at the time using *BayesTyper*. To run more samples, please execute *BayesTyper* in batches as described in the [batching wiki](https://github.com/bioinformatics-centre/BayesTyper/wiki/Executing-BayesTyper-on-sample-batches). Batching is currently not supported by the *snakemake* workflow - please let us know if you require this feature by filing a [feature request](https://github.com/bioinformatics-centre/BayesTyper/issues).

**Important:** Please note that Bayestyper currently does **not** support bgzip compressed vcf files, but only uncompressed and gzip compressed files.

### 1. Generation of variant candidates ###
Starting from a set of indexed, aligned reads (obtained e.g. using *BWA-MEM*):
1. For each sample, run [*HaplotypeCaller*](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php) to get standard mapping-based candidates
    * Marking duplicates, running Base Quality Recalibration or doing joint genotyping is **not** required
    * Faster alternative: [*Freebayes*](https://github.com/ekg/freebayes)
3. For each sample, run [*Platypus*](http://www.well.ox.ac.uk/platypus) to identify small and medium sized variants
4. For each sample, run [*Manta*](https://github.com/Illumina/manta) to identify candidates by *de novo* local assembly (important for detecting larger deletions and insertions). Convert allele IDs (e.g. \<DEL\>) in the Manta output to sequences using `bayesTyperTools convertAllele`
5. Combine variants across *all* samples, callers and the *variation prior* using `bayesTyperTools combine -o <candiate_variants> GATK:<gatk_sample1>.vcf,GATK:<gatk_sample2>.vcf,PLATYPUS:<platypus_sample1>.vcf,PLATYPUS:<platypus_sample2>.vcf,MANTA:<manta_sample1>.vcf,...,prior:<prior>.vcf`
   * **Note:** The source tag before the colon (e.g. GATK) only serves to identify the origin of the calls in the final callset - it can take any value.

### 2. Genotyping based on variant candidates ###
1. Count k-mers
   1. Run [KMC3](https://github.com/refresh-bio/KMC) on each sample (include singleton k-mers using `-ci1`)
      * **Note:** Default is fq(.gz) input - bam input is enabled using `-fbam`.
   2. Create a read k-mer bloom filter for each sample from the KMC3 output using `bayesTyperTools makeBloom -k <kmc_output_prefix> -p <num_threads>`
      * **Important:** The resulting bloom filter (*<sample_id>.bloom*) and the KMC3 output (*<sample_id>.kmc_pre* and *<sample_id>.kmc_suf*) must reside in the same directory and have the same prefix
2. Identify variant clusters: `bayesTyper cluster -v <candiate_variants>.vcf -s <samples>.tsv -g <ref_build>_canon.fa -d <ref_build>_decoy.fa -p <num_threads> -o <output_prefix>`
      * **Note:** This partitions the candidate variants into units written to separate directories (<output_prefix>_unit_1, <output_prefix>_unit_2 ...) - these can be processed independently e.g. on a cluster (supported by the snakemake workflow) followed by a simple concatenation operation using `bcftools concat` (see below).
      * **Important:** The `<samples>.tsv` file should contain one sample per row with columns \<sample_id\>, \<sex\> and \<kmc_output_prefix\> and no header ([example](http://people.binf.ku.dk/~lassemaretty/bayesTyper/bt_samples_example.tsv)).
      * **Important:** Data that is common to all units and necessary for genotyping is written to the <output_prefix>_cluster_data directory.
      * **Note:** Reference files (canon/decoy) are provided in the BayesTyper data bundle (see [Installation](https://github.com/bioinformatics-centre/BayesTyper#installation)).
3. Genotype variant clusters: `bayesTyper genotype -v <output_prefix>_unit_<unit_id>/variant_clusters.bin -c <output_prefix>_cluster_data -s <samples>.tsv -g <ref_build>_canon.fa -d <ref_build>_decoy.fa -o <output_prefix>_unit_<unit_id>/<output_prefix> -z -p <num_threads>`
      * **Note:** The genotype command also applies the default BayesTyper hard-filters by setting the variant FILTER status and the sample allele filter (SAF) format attribute. Please refer to the [filter wiki](https://github.com/bioinformatics-centre/BayesTyper/wiki/Filtering) for information about the filters used, how to changes the defaults up front and how to refilter the data after running the genotyping step.
      * **Important:** The filtering procedure **only** filters the genotypes ("./.") and set the genotype and variant level filter attributes. Hence, downstream tools **should prefilter on the variant quality and filter==PASS** (e.g. using *bcftools*) to obtain the filtered calls.
4. Concatenate units: `bcftools concat -O z -o <output_prefix>.vcf.gz <output_prefix>_unit_1/<output_prefix>.vcf.gz <output_prefix>_unit_2/<output_prefix>.vcf.gz ...`
    * **Important:** Unit arguments to bcftools should be in ascending order (unit_1, unit_2 ...) for the output to be properly sorted

## Computational requirements ##

|Number of samples|Coverage|Number of variant alleles|Max allele length (nts)|Number of threads|Wall time (h, single core)\*|Wall time (h, cluster)\**|Max memory (GB)|
|-|---|-----|-------|--|-|-|--|
|3|13x|21.4M|500,000|28|7|3|40|
|3|13x|64.4M|500,000|28|18|5|42|
|10|50x|11.7M|10,000|28|32|15|65|
|10|50x|61.1M|10,000|28|92|16|62|

\*The time estimates are for running `bayesTyper cluster` and `bayesTyper genotype` only. Expect <1h combined run-time per sample for counting k-mers by *KMC* and bloom filter creation by *bayesTyperTools*. All runs were done on a 64-bit Intel Xeon 2.00 GHz machine with 128 GB of memory.

\** `bayesTyper genotype` can be distributed across nodes on a cluster - between 2 and 11 nodes were used in this benchmark.

## Citing BayesTyper ##
The BayesTyper manuscript has been accepted for publication. Citation information will be updated on Monday June 18.

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
