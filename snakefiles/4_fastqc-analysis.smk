rule multiqc:
    input:
        expand(multiqc_raw + "/fastqc_report_htmls_zips/{stem}_R{R}_001_fastqc.zip", stem=stems, R=Rlist)
    output:
        multiqc_dir + "/fastqc_report_raw_reads.html",
        multiqc_dir + "/fastqc_report_raw_reads_data.zip"
    benchmark:
        bench + "/multi_qc.log"
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc {input} -o " + multiqc_dir + " -n fastqc_report_raw_reads -z"

rule fastqc:
    input:
        tmp + "/{stem}.fastq.gz"
    output:
        multiqc_raw + "/fastqc_report_htmls_zips/{stem}_fastqc.zip",
        multiqc_raw + "/fastqc_report_htmls_zips/{stem}_fastqc.html"
    benchmark:
        bench + "/fastqc_{stem}.log"
    params:
        kmer_size = config["FASTQC__kmer_size"],
        out_dir = multiqc_raw + "/fastqc_report_htmls_zips/"
    conda:
        "../envs/fastqc.yaml"
    threads:
        config_threads["FASTQC"]
    shell:
        "fastqc {input} -o {params.out_dir} -t {threads} -k {params.kmer_size}"

rule multiqc_processed:
    input:
        get_multiqc_pre
    output:
        multiqc_dir + "/fastqc_report_processed_reads.html",
        multiqc_dir + "/fastqc_report_processed_reads_data.zip"
    benchmark:
        bench + "/multi_qc_pre.log"
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc {input} -o " + multiqc_dir + " -n fastqc_report_processed_reads --force -z"
