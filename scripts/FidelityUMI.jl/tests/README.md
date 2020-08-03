# Getting test `BAM` files

* Compress test fasta files:
```bash
    gzip -k scripts/FidelityUMI.jl/tests/1_duplicated*
    gzip -k scripts/FidelityUMI.jl/tests/2_duplicated*
```

* To generate tes `BAM` files run `ngsanalysis` pipeline:  
```bash
   snakemake --use-conda --configfile scripts/Fidelity.jl/tests/umi_fidelity.yaml umi_fidelity --nt
   snakemake --use-conda --configfile scripts/Fidelity.jl/tests/umi_fidelity2.yaml umi_fidelity --nt
```

* Created `BAM` files needed can be found in temporary directories:
```python
   tmp + "/{stem}_umi_sort_by_name.bam"
```

* Exact list of files:
   1. 1_duplicated_umi_sort_by_name.bam
   2. 2_duplicated_umi_sort_by_name.bam

* These files can be copied to test directory if needed.

# Generating test fastq files with `dwgsim`

1. Install your environment `dwgsim`. In this case `conda`:
```bash
    conda install -c bioconda dwgsim
```

2. Generate pair end reads with known mutations where reference is `references/fidelity/Lambda_chr.fasta` and output is `dwgsim_umi_fidelity_errors`:
```bash
    dwgsim -N 20 -r 0.02  {reference} {output}
```
This will create 20 reads from provided reference with error frequency ~ 0.02.

3. To make completely overlapping pairs, first rename reads in R2 file by replacing part of the read name which exists in all read names, for example ID:
```bash
    sed -i 's/@Enterobacteria_phage_lambda/@Enterobacteria_phage_lambda_read2/g' dwgsim_umi_fidelity_errors_renamed.bwa.read2.fastq
```

4. Compress both `fastq` files:
```bash
    gzip -k {read1_file} {read2_file}
```

5. Get reverse complement for both `fastq` files:
```bash
     julia scripts/fastq_reverse_complement.jl -f {read1_file.fastq.gz} -o {read1_reverse_complement_file.fastq.gz}
     julia scripts/fastq_reverse_complement.jl -f {read2_file.fastq.gz} -o {read2_reverse_complement_file.fastq.gz}
```
At this point you have 4 `fastq.gz` files

6. Merge read1 with read2 files:
```bash
    cat {read1_file.fastq.gz} {read2_file.fastq.gz} > your_filename_L001_R1_001.fastq.gz
    cat {read1_reverse_complement_file.fastq.gz} {read2_reverse_complement_file.fastq.gz} > your_filename_L001_R2_001.fastq.gz
```

7. Calculate error rate with `paired` algorithm:
```bash
   snakemake --use-conda --configfile scripts/Fidelity.jl/tests/umi_fidelity_lambda.yaml umi_fidelity
```

8. Add more background errors:
```bash
    dwgsim -r 0.001 -H -S 2 {reference} {output2}
    gzip -k {output2}.bwa.read1.fastq {output2}.bwa.read2.fastq
    cat {output2}.bwa.read1.fastq.gz >> {output2}_L001_R1_001.fastq.gz
    cat {output2}.bwa.read2.fastq.gz >> {output2}_L001_R2_001.fastq.gz
```

9. Calculate error rate with `paired` algorithm:
```bash
   snakemake --use-conda --configfile scripts/Fidelity.jl/tests/umi_fidelity_lambda.yaml umi_fidelity
```
