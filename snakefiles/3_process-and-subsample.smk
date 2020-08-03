rule umi_tools_extract:
    input:
        rules.get_fastq1.output,
        rules.get_fastq2.output,
    output:
        temp(tmp + "/{stem}_R1_001_umi.fastq.gz"),
        temp(tmp + "/{stem}_R2_001_umi.fastq.gz"),
        logs + "/EXTRACT-UMI/{stem}_umi_extract.log"
    params:
        umi1 = config["UMI_TOOLS__umi_structure"],
        add = config["UMI_TOOLS__additional_extract_params"]
    conda:
        "../envs/umi-tools.yaml"
    benchmark:
        bench + "/extract_umi_{stem}.log"
    shell:
        "umi_tools extract --stdin={input[0]} --bc-pattern={params.umi1} " +
        "--log={output[2]} --stdout={output[0]} " +
        "--read2-in={input[1]} --read2-out={output[1]} " +
        "{params.add}"

rule bbduk:
    input:
        get_trimm_samples,
        adapters_bbduk
    output:
        temp(tmp + "/{stem}_R1_001_bbduk.fastq.gz"),
        temp(tmp + "/{stem}_R2_001_bbduk.fastq.gz"),
        temp(tmp + "/{stem}_trim.settings")
    benchmark:
        bench + "/trim_adapters_{stem}.log"
    params:
        ktrim = config["TRIMMING__bbduk__ktrim"],
        k = config["TRIMMING__bbduk__k"],
        mink = config["TRIMMING__bbduk__mink"],
        hdist = config["TRIMMING__bbduk__hdist"],
        minlength = config["TRIMMING__bbduk__minlength"],
        maxns = config["TRIMMING__bbduk__maxns"],
        qtrim = config["TRIMMING__bbduk__qtrim"],
        trimq = config["TRIMMING__bbduk__trimq"],
        additional_params = config["TRIMMING__bbduk__additional_params"],
        mem = mem_java
    conda:
        "../envs/bbmap.yaml"
    threads:
        config_threads["BBDUK"]
    shell:
        "bbduk.sh in={input[0]} in2={input[1]} out={output[0]} out2={output[1]} " +
        "ref={input[2]} threads={threads} ordered=t -Xmx{params.mem}g " +
        "ktrim={params.ktrim} k={params.k} mink={params.mink} hdist={params.hdist} " +
        "minlength={params.minlength} maxns={params.maxns} {params.additional_params} " +
        "qtrim={params.qtrim} trimq={params.trimq} " +
        "2>{output[2]}; "

rule adapter_removal:
    input:
        get_trimm_samples,
        adapters_adapt
    output:
        temp(tmp + "/{stem}_R1_001_adapt.fastq.gz"),
        temp(tmp + "/{stem}_R2_001_adapt.fastq.gz"),
        temp(tmp + "/{stem}_trim.settings")
    benchmark:
        bench + "/trim_adapters_{stem}.log"
    params:
        quality_limit=config["TRIMMING__AdapterRemoval__quality_limit"],
        maximum_n_count_in_reads=config["TRIMMING__AdapterRemoval__maximum_n_count_in_reads"],
        max_errors=config["TRIMMING__AdapterRemoval__max_errors"],
        min_length=config["TRIMMING__AdapterRemoval__min_read_length"]
    conda:
        "../envs/adapterremoval.yaml"
    threads:
        config_threads["ADAPTERREMOVAL"]
    shell:
        "AdapterRemoval --threads {threads} " +
        " --maxns {params.maximum_n_count_in_reads} " +
        " --minlength {params.min_length} --trimqualities " +
        " --minquality {params.quality_limit} --adapter-list {input[2]} " +
        " --output1 {output[0]} --output2 {output[1]} " +
        " --file1 {input[0]}  --file2 {input[1]} --settings {output[2]} " +
        " --outputcollapsed " + tmp + "/{wildcards.stem}.collapsed" +
        " --outputcollapsedtruncated " + tmp + "/{wildcards.stem}.collapsed.truncated" +
        " --singleton " + tmp + "/{wildcards.stem}_singleton.truncated.fastq.gz" +
        " --discarded " + tmp + "/{wildcards.stem}_discarded.fastq.gz " +
        " --gzip --barcode-mm 2 --barcode-mm-r1 1 --barcode-mm-r2 1"


rule plot_clean_bases_bbduk:
    input:
        expand(tmp + "/{stem}_trim.settings",stem=stems)
    output:
        trim_dir + "/clean_bases_bbduk.png",
        trim_dir + "/clean_bases_bbduk.csv"
    shell:
        "python "+CURRENT_PATH+"/scripts/plot_clean_bases_bbduk.py --input {input} --out {output[0]}"

rule make_adapters_for_bbduk:
    input:
        config["TRIMMING__adapters"]
    output:
        tmp + "/bbduk_adapters.fa"
    shell:
        "awk 'BEGIN{{a=1;OFS=\"\\n\"}}{{print \">Adapter_\"a,$1;a++;if($2){{print \">Adapter_\"a,$2;a++}}}}' " +
        "{input} > {output}"



rule bbduk_force_trim:
    input:
        get_force_trimm_samples
    output:
        temp(tmp + "/{stem}_R1_001_2.fastq.gz"),
        temp(tmp + "/{stem}_R2_001_2.fastq.gz"),
        temp(tmp + "/{stem}_trim_force.settings")
    benchmark:
        bench + "/trim_adapters2_{stem}.log"
    params:
        minlength = config["TRIMMING__bbduk__minlength"],
        force_options = config["TRIMMING__bbduk__FORCE_TRIM_OPTIONS"],
        mem = mem_java
    conda:
        "../envs/bbmap.yaml"
    threads:
        config_threads["BBDUK"]
    shell:
        "bbduk.sh in={input[0]} in2={input[1]} out={output[0]} out2={output[1]} " +
        "threads={threads} ordered=t -Xmx{params.mem}g " +
        "minlength={params.minlength} {params.force_options} " +
        "2>{output[2]}; "

rule plot_clean_bases_bbduk2:
    input:
        expand(tmp + "/{stem}_trim_force.settings", stem=stems)
    output:
        trim_dir + "/force_trim_bases_bbduk.png",
        trim_dir + "/force_trim_bases_bbduk.csv"
    shell:
        "python "+CURRENT_PATH+"/scripts/plot_clean_bases_bbduk.py --input {input} --out {output[0]}"

