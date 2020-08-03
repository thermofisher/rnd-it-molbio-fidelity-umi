rule samtools_flagstat1:
    input:
        aligned_dir + "/{stem}_subSort.bam",
        aligned_dir + "/{stem}_subSort.bam.bai"
    output:
        out + "/FIDELITY/{stem}_subSort_flagstat.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools flagstat {input[0]} > {output}"

rule samtools_fix_mates:
    input:
        aligned_dir + "/{stem}_subSort.bam"
    output:
        temp(tmp + "/{stem}_fixed_mates.bam")
    conda:
        "../envs/samtools.yaml"
    threads:
        config_threads["SAMTOOLS_FIX_MATES"]
    shell:
        "samtools view {input} -o {output} -b -F 4 -F 8 -F 256 -F 2048"

rule umi_tools_group_reads:
    input:
        tmp + "/{stem}_fixed_mates.bam",
        tmp + "/{stem}_fixed_mates.bam.bai"
    output:
        temp(tmp + "/{stem}_grouped.bam")
    log:
        logs + "/UMI-GROUP/{stem}.log"
    params:
        haming = config["UMI_TOOLS__group_haming_distance"]
    conda:
        "../envs/umi-tools.yaml"
    threads:
        config_threads["UMI_GROUP"]
    shell:
        "umi_tools umi_tools group -I {input[0]} -S {output} -L {log} " +
        "--paired --method adjacency --output-bam " +
        "--edit-distance-threshold {params.haming}"

rule samtools_flagstat2:
    input:
        calculate_fidelity_sort_by_name_input
    output:
        out + "/FIDELITY/{stem}_flagstat_sort_by_name.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools flagstat {input[0]} > {output}"

rule sort_by_name:
    input:
        calculate_fidelity_sort_by_name_input
    output:
        temp(tmp + "/{stem}_sort_by_name.bam")
    params:
        temp_dir_sambamba = tmp + '/SAMBAMBA'
    conda:
        "../envs/sambamba.yaml"
    threads:
        config_threads["SAMBAMBA_sort"]
    shell:
        "sambamba sort -n --tmpdir={params.temp_dir_sambamba} -p " +
        "-t{threads} -o {output} {input}"

rule get_fidelity:
    input:
        tmp + "/{stem}_sort_by_name.bam"
    output:
        umi_fidelity_output
    log:
        logs + "/CALCULATE-FIDELITY-UMI-BASED/{stem}.log"
    benchmark:
        bench + "/calculate_fidelity_{stem}.txt"
    params:
        p = umi_fidelity_params
    threads:
        config_threads["FIDELITY_UMI"]
    shell:
        "julia "+CURRENT_PATH+"/scripts/FidelityUMI.jl/scripts/FidelityUmi.jl " +
        "-b {input} -o {output[0]} {params.p} &> {log}"

rule concatenate_results:
    input:
        expand(tmp + "/{stem}_fidelity.csv", stem=stems)
    output:
        out + "/FIDELITY/fidelity.csv"
    shell:
        "julia "+CURRENT_PATH+"/scripts/join_csv.jl {input} {output}"

rule merge_fidelity_transitversions:
    input:
        expand(tmp + "/{stem}_fidelity_transitversions.csv", stem=stems)
    output:
        out + "/FIDELITY/transitversions.csv"
    shell:
        "julia "+CURRENT_PATH+"/scripts/join_csv.jl {input} {output}"

rule merge_fidelity_transvers_types:
    input:
        expand(tmp + "/{stem}_fidelity_transvers_types.csv", stem=stems)
    output:
        out + "/FIDELITY/transvers_types.csv"
    shell:
        "julia "+CURRENT_PATH+"/scripts/join_csv.jl {input} {output}"

rule merge_fidelity_mutations:
    input:
        expand(tmp + "/{stem}_fidelity_passing_mutations.csv", stem=stems)
    output:
        out + "/FIDELITY/passing_mutations.csv"
    shell:
        "julia "+CURRENT_PATH+"/scripts/join_csv.jl {input} {output}"

rule merge_fidelity_single_mutations:
    input:
        expand(tmp + "/{stem}_fidelity_passing_single_mutations.csv", stem=stems)
    output:
        out + "/FIDELITY/passing_single_mutations.csv"
    shell:
        "julia "+CURRENT_PATH+"/scripts/join_csv.jl {input} {output}"

rule merge_fidelity_translation:
    input:
        expand(tmp + "/{stem}_fidelity_passing_aa_mutations.csv", stem=stems)
    output:
        out + "/FIDELITY/aa_mutations.csv"
    shell:
        "julia "+CURRENT_PATH+"/scripts/join_csv.jl {input} {output}"

rule merge_fidelity_single_translation:
    input:
        expand(tmp + "/{stem}_fidelity_passing_single_aa_mutations.csv", stem=stems)
    output:
        out + "/FIDELITY/single_aa_mutations.csv"
    shell:
        "julia "+CURRENT_PATH+"/scripts/join_csv.jl {input} {output}"

rule merge_fidelity_frameshifts:
    input:
        expand(tmp + "/{stem}_fidelity_passing_frameshifts.csv", stem=stems)
    output:
        out + "/FIDELITY/frameshifts.csv"
    shell:
        "julia "+CURRENT_PATH+"/scripts/join_csv.jl {input} {output}"

rule plot_frequencies_fidelity_based:
    input:
        tmp + "/{stem}_fidelity_distributions.csv"
    output:
        out + "/FIDELITY/{stem}_fidelity_errors_distributions.png"
    shell:
        "python "+CURRENT_PATH+"/scripts/plot_hist.py --input {input} " +
        "--output {output} --log --hist-type stepfilled --bins 100 " +
        "--y-label Count --x-label Error_frequency"

rule plot_umi_size_distributions:
    input:
        tmp + "/{stem}_fidelity_size_distribution.csv"
    output:
        out + "/FIDELITY/{stem}_fidelity_size_distribution.png"
    shell:
        "python "+CURRENT_PATH+"/scripts/plot_hist.py --input {input} " +
        "--output {output} --log --hist-type stepfilled --bins 100 "+
        "--y-label Count --x-label Cluster_size "

rule plot_umi_error_distributions_read:
    input:
        tmp + "/{stem}_fidelity_error_distributions_per_read.csv"
    output:
        out + "/FIDELITY/{stem}_fidelity_errors_per_read_distribution.png"
    shell:
        "python "+CURRENT_PATH+"/scripts/plot_hist.py --input {input} " +
        "--output {output} --hist-type stepfilled --bins 100 --use-columns-as-array " +
        "--y-label Frequency --x-label Position,bp "

rule plot_cumulative_err:
    input:
        out + "/FIDELITY/passing_single_mutations.csv"
    output:
        out + "/FIDELITY/cumulative.png",
        out + "/FIDELITY/cumulative_graph_data.csv"
    shell:
        "python "+CURRENT_PATH+"/scripts/plot_cumulative.py --input {input} " +
        "--output {output[0]} --x-name ClustOfVariants --y-name Count " +
        "--z-name ConfidentCoverage --sample-name Sample --chromosome-name Chromosome " +
        "--y-label 'Fraction of confident bases' --x-label 'Frequency of errors'"

rule plot_scatter_errors:
    input:
        out + "/FIDELITY/passing_single_mutations.csv"
    output:
        out + "/FIDELITY/single_mutations.png",
        out + "/FIDELITY/single_mutations_graph_data.csv",
    shell:
        "python "+CURRENT_PATH+"/scripts/plot_errors_per_chromosome.py --input {input} " +
        "--output {output[0]} --x-name ClustOfVariants --y-name Count " +
        "--z-name ConfidentCoverage --sample-name Sample --chromosome-name Chromosome " +
        "--matplotlib-kwargs alpha=0.25 marker=o " +
        "--y-label 'Error frequency, errors/reads' --x-label 'position, bp'"
