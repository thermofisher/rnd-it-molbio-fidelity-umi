rule get_files_to_cat:
    output:
        tmp + "/{stem}R1_input_list.sh",
        tmp + "/{stem}R2_input_list.sh"
    run:
        get_sample_files_script(wildcards.stem)

rule get_fastq1:
    input:
        tmp + "/{stem}R1_input_list.sh"
    output:
        tmp + "/{stem}_R1_001.fastq.gz"
    threads:
        config_threads["get_fastq"]
    shell:
        "./{input} > {output}"

rule get_fastq2:
    input:
        tmp + "/{stem}R2_input_list.sh"
    output:
        tmp + "/{stem}_R2_001.fastq.gz"
    threads:
        config_threads["get_fastq"]
    shell:
        "./{input} > {output}"
