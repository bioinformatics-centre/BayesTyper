########################################
# STAGE 2: GENOTYPE CANDIDATE VARIANTS #
########################################

rule kmc_count_kmers:
    input:
        wildscardsToReadmap
    output:
        kmc_pre="kmer_counts/{sample_id}.kmc_pre",
        kmc_suf="kmer_counts/{sample_id}.kmc_suf",
        tmp="kmer_counts/tmp/{sample_id}/"
    log:
        "kmer_counts/{sample_id}_kmc.log"
    params:
        out_prefix="kmer_counts/{sample_id}",
        runtime="4:00:00"
    shell:
        "{config[kmc_path]} -k55 -ci1 -fbam {input} {params.out_prefix} {output.tmp} > {log} 2>&1"

rule bayestyper_make_bloom:
    input:
        "kmer_counts/{sample_id}.kmc_pre",
        "kmer_counts/{sample_id}.kmc_suf"
    output:
        "kmer_counts/{sample_id}.bloomData",
        "kmer_counts/{sample_id}.bloomMeta"
    log:
        "kmer_counts/{sample_id}_bloom.log"
    params:
        out_prefix="kmer_counts/{sample_id}",
        runtime="4:00:00"
    shell:
        "{bayestyper_tools_path} makeBloom -k {params.out_prefix} > {log} 2>&1"

rule bayestyper_make_samples_file:
    output:
        "bayestyper_samples.tsv"
    params:
        runtime="00:05:00"
    run:
        with open(output[0], "w") as bayestyper_samples_file:
            for sample_idx in range(samples.shape[0]):
                bayestyper_samples_file.write(samples.iloc[sample_idx,:].sample_id + "\t" + samples.iloc[sample_idx,:].sex + "\t" + "kmer_counts/" + samples.iloc[sample_idx,:].sample_id + "\n")

rule bayestyper_cluster:
    input:
        expand("kmer_counts/{sample_id}.bloomData", sample_id=samples.loc[:,"sample_id"]),
        expand("kmer_counts/{sample_id}.bloomMeta", sample_id=samples.loc[:,"sample_id"]),
        variants="candidates/candidates.vcf.gz",
        samples="bayestyper_samples.tsv"
    output:
        units=dynamic("bayestyper/bayestyper_unit_{unit_id}/variant_clusters.bin")
    params:
        out_prefix ="bayestyper/bayestyper",
        runtime="12:00:00",
        log="bayestyper/cluster.log"
    threads: 16
    shell:
        "{bayestyper_path} cluster -v {input.variants} -s {input.samples} -g {bayestyper_reference_canon} -d {bayestyper_reference_decoy} -p {threads} -o {params.out_prefix} > {params.log} 2>&1"

rule bayestyper_genotype:
    input:
        kmc_pre=expand("kmer_counts/{sample_id}.kmc_pre", sample_id=samples.loc[:,"sample_id"]),
        kmc_suf=expand("kmer_counts/{sample_id}.kmc_suf", sample_id=samples.loc[:,"sample_id"]),
        samples="bayestyper_samples.tsv",
        unit="bayestyper/bayestyper_unit_{unit_id}/variant_clusters.bin"
    output:
        genotypes="bayestyper/bayestyper_unit_{unit_id}/bayestyper.vcf.gz",
        kmer_coverage_file="bayestyper/bayestyper_unit_{unit_id}/bayestyper_genomic_parameters.txt"
    log:
        "bayestyper/bayestyper_unit_{unit_id}/genotype.log"
    params:
        cluster_data_dir="bayestyper/bayestyper_cluster_data/",
        out_prefix="bayestyper/bayestyper_unit_{unit_id}/bayestyper",
        runtime="12:00:00"
    threads: 16
    shell:
        "{bayestyper_path} genotype -v {input.unit} -s {input.samples} -c {params.cluster_data_dir} -g {bayestyper_reference_canon} -d {bayestyper_reference_decoy} -p {threads} -z -o {params.out_prefix} > {log} 2>&1"

rule bcftools_concat_units:
    input:
        dynamic("bayestyper/bayestyper_unit_{unit_id}/bayestyper.vcf.gz")
    output:
        "bayestyper/bayestyper.vcf.gz"
    params:
        runtime="12:00:00"
    run:
        sorted_input = sorted(input)
        print(sorted_input)
        shell("{config[bcftools_path]} concat -O z -o {output} {sorted_input}")
