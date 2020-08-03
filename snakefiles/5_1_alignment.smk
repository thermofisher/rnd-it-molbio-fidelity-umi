rule bwa_mem:
    input:
        get_processed_samples
    output:
        temp(tmp + "/{stem}_aligned.bam"),
        temp(tmp + "/{stem}_alignmnt.log")
    benchmark:
        bench + "/bwa_{stem}.log"
    params:
        reference = config["MAIN__REFERENCE"],
        non_default_bwa_params = config["BWA__non_default_params"]
    conda:
        "../envs/bwa.yaml"
    threads:
        config_threads["BWA"]
    shell:
        "bwa mem -t {threads} {params.reference} "+
        "{params.non_default_bwa_params} {input} "+
        "| samtools view -@ 3 -bS - > {output[0]}  2> {output[1]}"

rule bwa_index:
    input:
        config["MAIN__REFERENCE"]
    output:
        config["MAIN__REFERENCE"] + ".bwt"
    benchmark:
        bench + "/bwa_index.log"
    conda:
        "../envs/bwa.yaml"
    shell:
        "bwa index -a bwtsw {input}"

rule sambamba_sort:
    input:
        tmp + "/{stem}_aligned.bam"
    output:
        aligned_dir + "/{stem}_subSort.bam"
    benchmark:
        bench + "/sortMaped_{stem}.log"
    params:
        temp_dir_sambamba = tmp + '/SAMBAMBA'
    conda:
        "../envs/sambamba.yaml"
    threads:
        config_threads["SAMBAMBA_sort"]
    shell:
        "sambamba sort --tmpdir={params.temp_dir_sambamba} -p " +
        "-t{threads}  -o {output[0]} {input[0]}"

rule sambamba_index:
    input:
        tmp + "/{stem}.bam"
    output:
        temp(tmp + "/{stem}.bam.bai")
    benchmark:
        bench + "/sambamba_index_{stem}.log"
    threads:
        config_threads["SAMBAMBA_index"]
    conda:
        "../envs/sambamba.yaml"
    shell:
        "sambamba index -p -t{threads} {input[0]}"
