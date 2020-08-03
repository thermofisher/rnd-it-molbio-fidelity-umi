import collections
import logging
import logging.config
import os
import time
import yaml
import glob
import re
import subprocess
import sys
import getpass

# ------------------------------ Include Snakefiles --------------------------- #

include: "./snakefiles/0_1_utilities.smk"
include: "./snakefiles/0_2_configurations.smk"
include: "./snakefiles/1_download.smk"
include: "./snakefiles/3_process-and-subsample.smk"
include: "./snakefiles/4_fastqc-analysis.smk"
include: "./snakefiles/5_1_alignment.smk"
include: "./snakefiles/7_1_DNA-alignment-analysis.smk"

# ------------------------------ Target description --------------------------- #

target = [
    logs + "/GIT-LOG/commit_used.log",
    logs + "/config.yaml",
    multiqc_dir + "/fastqc_report_raw_reads.html",
    expand(multiqc_raw + "/fastqc_report_htmls_zips/{stem}_R{R}_001_fastqc.zip", stem=stems, R=Rlist),
    multiqc_dir + "/fastqc_report_raw_reads_data.zip",
    multiqc_dir + "/fastqc_report_processed_reads.html",
    expand(multiqc_raw + "/fastqc_report_htmls_zips/{stem}_R{R}_001_2_fastqc.zip", stem=stems, R=Rlist),
    multiqc_dir + "/fastqc_report_processed_reads_data.zip",
    trim_dir + "/clean_bases_bbduk.png",
    trim_dir + "/clean_bases_bbduk.csv"
]


umi_fidelity_rule_targets = [
    expand(out + "/FIDELITY/{stem}_fidelity_errors_distributions.png", stem=stems),
    expand(out + "/FIDELITY/{stem}_fidelity_size_distribution.png", stem=stems),
    expand(out + "/FIDELITY/{stem}_flagstat_sort_by_name.txt", stem=stems),
    expand(out + "/FIDELITY/{stem}_subSort_flagstat.txt", stem=stems),
    expand(out + "/FIDELITY/{stem}_fidelity_errors_per_read_distribution.png", stem=stems),
    out + "/FIDELITY/fidelity.csv",
    out + "/FIDELITY/transitversions.csv",
    out + "/FIDELITY/transvers_types.csv",
    out + "/FIDELITY/cumulative.png",
    out + "/FIDELITY/single_mutations.png",
    out + "/FIDELITY/cumulative_graph_data.csv",
    out + "/FIDELITY/single_mutations_graph_data.csv",
]

if config["FIDELITY__report_mutations"]:
    umi_fidelity_rule_targets += [
        out + "/FIDELITY/passing_mutations.csv",
        out + "/FIDELITY/passing_single_mutations.csv"
    ]
if config["FIDELITY__genes_csv"]:
    umi_fidelity_rule_targets += [
        out + "/FIDELITY/aa_mutations.csv",
        out + "/FIDELITY/frameshifts.csv",
        out + "/FIDELITY/single_aa_mutations.csv"
    ]
target.extend(umi_fidelity_rule_targets)

rule all:
    input: target

rule umi_fidelity:
    input:
        umi_fidelity_rule_targets

rule trim:
    input:
        expand(tmp + "/{stem}_R1_001_2.fastq.gz", stem=stems),
        expand(tmp + "/{stem}_R2_001_2.fastq.gz", stem=stems),
        trim_dir + "/clean_bases_bbduk.csv",
        trim_dir + "/force_trim_bases_bbduk.csv",
        trim_dir + "/clean_bases_bbduk.png",
        trim_dir + "/force_trim_bases_bbduk.png",

rule multiqc_reads:
    input:
        multiqc_dir + "/fastqc_report_processed_reads.html",
        multiqc_dir + "/fastqc_report_raw_reads.html"

rule get_index:
    input:
        index_file

rule download:
    input:
        expand(tmp + "/{stem}_R1_001.fastq.gz", stem=stems),
        expand(tmp + "/{stem}_R2_001.fastq.gz", stem=stems)

rule git_log:
    output:
        logs + "/GIT-LOG/commit_used.log"
    shell:
        "echo $(git log) |  cut -d ' ' -f 2 &>> {output}"

rule write_config:
    output:
        logs + "/config.yaml"
    run:
        with open(logs + "/config.yaml", 'w') as outfile:
            yaml.dump(config, outfile, default_flow_style=False)