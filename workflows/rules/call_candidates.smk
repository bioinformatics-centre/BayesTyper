#########################################
# STAGE 1: ENUMERATE CANDIDATE VARIANTS #
#########################################
rule gatk_call_variants:
    input:
        wildscardsToReadmap
    output:
        temp("gatk/{sample_id}.vcf")
    log:
        "gatk/{sample_id}.log"
    params:
        runtime = "12:00:00"
    threads: 16
    shell:
        "{config[java_path]} -jar {config[gatk_path]} -T HaplotypeCaller -nct {threads} -R {reference} -I {input} --genotyping_mode DISCOVERY -stand_call_conf 10 -o {output} > {log} 2>&1"

rule platypus_call_variants:
    input:
        wildscardsToReadmap
    output:
        temp("platypus/{sample_id}.vcf")
    log:
        "platypus/{sample_id}.log"
    params:
        runtime = "12:00:00"
    threads: 16
    shell:
        "{config[python2.7_path]} {config[platypus_path]} callVariants --bamFiles={input} --refFile={reference} --output={output} --logFileName={log} --nCPU {threads} --assemble=1 --assembleBrokenPairs=1"

rule platypus_make_contigs: # Platypus vcf output lacks contigs in the header - these are required downstream
    output:
        temp("platypus/contigs.txt")
    log:
        "platypus/makeContigs.log"
    params:
        runtime = "00:15:00"
    run:
        import sys
        contigs_file = open(output[0], "w")
        with open(reference_dict) as ref_dict_file:

            for line in ref_dict_file:

                if line[:3] == "@SQ":
                    contigs_file.write("##contig=<ID=" + line.split("\t")[1].split(":")[1] + ",length=" + line.split("\t")[2].split(":")[1] + ">\n")

        contigs_file.close()

rule platypus_norm:
    input:
        variants="platypus/{sample_id}.vcf",
        contigs="platypus/contigs.txt"
    output:
        temp("platypus/{sample_id}_norm.vcf.gz"),
    params:
        runtime = "4:00:00"
    shell:
        "sed -e '/##platypusOptions/r {input.contigs}' {input.variants} | {config[bcftools_path]} norm -c x -f {reference} -m -any -w 1000000 -O v | gzip -c > {output}" # sed command required to fix missing contigs in Platypus vcf

rule manta_prepare_workflow:
    input:
        wildscardsToReadmap
    output:
        dir="manta/{sample_id}",
        file="manta/{sample_id}/runWorkflow.py"
    log:
        "manta/{sample_id}_prepare.log"
    params:
        runtime = "4:00:00"
    shell:
        "{config[python2.7_path]} {config[manta_config_path]} --bam {input} --referenceFasta {reference} --runDir {output.dir} > {log} 2>&1" # Snakemake requires python3 - but manta requires python2 - therefore explicit call to python2

rule manta_run_workflow:
    input:
        "manta/{sample_id}/runWorkflow.py"
    output:
        temp("manta/{sample_id}/results/variants/candidateSV.vcf.gz")
    log:
        "manta/{sample_id}_run.log"
    params:
        runtime = "4:00:00",
        gunzipped_output = "manta/{sample_id}/results/variants/candidateSV.vcf"
    shell:  # gunzip -> gzip roundtrip required as BayesTyper v1.3 cannot read bgzip compressed files due to issues with the Boost iostreams gzip library
        """
        {input} -m local -g 63 > {log}
        gunzip {output}
        gzip {params.gunzipped_output}
        """

rule manta_conv_all_id:
    input:
        "manta/{sample_id}/results/variants/candidateSV.vcf.gz"
    output:
        temp("manta/{sample_id}/results/variants/candidateSV_converted.vcf.gz")
    log:
        "manta/{sample_id}_convert.log"
    params:
        runtime = "8:00:00"
    run:
        output_prefix = output[0].split(".")[0]
        shell("{bayestyper_tools_path} convertAllele -v {input} -g {reference} -z -o {output_prefix} > {log}")

rule manta_norm:
    input:
        "manta/{sample_id}/results/variants/candidateSV_converted.vcf.gz"
    output:
        temp("manta/{sample_id}/results/variants/candidateSV_converted_norm.vcf.gz")
    log:
        "manta/{sample_id}_norm.log"
    params:
        runtime = "4:00:00"
    shell:
        "{config[bcftools_path]} norm -c x -f {reference} -m -any -w 1000000 -O v {input} | gzip -c > {output} 2> {log}"

rule bayestyper_combine_variants:
    input:
        gatk=expand("gatk/{sample_id}.vcf", sample_id=samples.loc[:,"sample_id"]),
        platypus=expand("platypus/{sample_id}_norm.vcf.gz", sample_id=samples.loc[:,"sample_id"]),
        manta=expand("manta/{sample_id}/results/variants/candidateSV_converted_norm.vcf.gz", sample_id=samples.loc[:,"sample_id"]),
        prior=bayestyper_variant_prior_path
    output:
        protected("candidates/candidates.vcf.gz")
    log:
        "candidates/combine.log"
    params:
        runtime = "4:00:00"
    run:
        input_vcf_str = ",".join(["gatk:" + x for x in input.gatk] + ["manta:" + x for x in input.manta] + ["prior:" + input.prior])
        output_prefix = output[0].split(".")[0]
        shell("{bayestyper_tools_path} combine -v {input_vcf_str} -z -o {output_prefix} > {log}")
