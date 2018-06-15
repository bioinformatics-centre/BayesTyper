import pandas, os, sys

configfile: "config.yaml"

def wildscardsToReadmap(wildcards):

    return samples.loc[wildcards.sample_id].read_map

def getPrefix(path, num_suffixes = 1):

    return ".".join(path.split(".")[:-num_suffixes])

samples = pandas.read_table(config["samples_file"]).set_index("sample_id", drop=False)
if not config["reference_build"] in ["GRCh37", "GRCh38"]:

    sys.exit("ERROR: Reference build must be GRCh37 or GRCh38!")

reference_prefix = os.path.join(config["bayestyper_data_bundle_dir"], config["reference_build"])
reference = reference_prefix + ".fa"
reference_dict = reference_prefix + ".dict"
bayestyper_reference_canon = reference_prefix + "_canon.fa"
bayestyper_reference_decoy = reference_prefix + "_decoy.fa"
bayestyper_variant_prior_path = os.path.join(config["bayestyper_data_bundle_dir"],"SNP_dbSNP150common_SV_1000g_dbSNP150all_GDK_GoNL_GTEx_" + config["reference_build"] + ".vcf.gz")

bayestyper_path = os.path.join(config["bayestyper_bin_dir"], "bayesTyper")
bayestyper_tools_path = os.path.join(config["bayestyper_bin_dir"], "bayesTyperTools")

rule all:
    input:
        "bayestyper/bayestyper.vcf.gz"

include: "rules/call_candidates.smk"
include: "rules/estimate_genotypes.smk"
