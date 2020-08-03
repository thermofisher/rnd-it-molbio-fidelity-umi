#!/usr/bin/env julia

using ArgParse.ArgParseSettings
using ArgParse.@add_arg_table
using ArgParse.parse_args
using DataFrames
using CSV
using StatsBase

NEW_PATH=join(split(Base.source_path(),"/")[1:end-3],"/")
push!(LOAD_PATH, NEW_PATH)
using FidelityUMI


function main(args)
    arg_parse_settings = ArgParseSettings(description="Calculates N of variants based on UMIs groups and filtering sequencing errors.
                                                      UMI_tools should be used for tagging readname and/or grouping.")
    @add_arg_table arg_parse_settings begin

        "--bam", "-b"
            arg_type = String
            help = "name of bam file."
            required = true
            dest_name = "bam"

        "--output", "-o"
            arg_type = String
            help = "Output csv file"
            default = "output.csv"
            dest_name = "out"

        "--reference", "-r"
            arg_type = String
            help = "Reference fasta."
            required = true
            dest_name = "ref"

        "--alignment-quality", "-q"
            arg_type = Int64
            help = "Minimum alignment quality of a read."
            default = 0
            dest_name = "alignQ"

        "--cluster-size", "-c"
            arg_type = Int64
            help = "How many read pairs in cluster should be considered."
            default = 0
            dest_name = "clust_size"

        "--mut-per-read", "-m"
            arg_type = Int64
            help = "Allowed max number of mutations in read."
            default = 3
            dest_name = "mut_n_read"

        "--consider-as-ref-error", "-e"
            arg_type = Float64
            help = "Consider mutations at position as reference error if fraction is >= x."
            default = 1.0
            dest_name = "ref_error"

        "--base-quality", "-Q"
            arg_type = Int64
            help = "Base quality [>=]."
            default = 30
            dest_name = "q_limit"

        "--umi-string", "-s"
            arg_type = String
            help = "Umi place holder string."
            default = "_"
            dest_name = "umi_s"

        "--umi-place", "-p"
            arg_type = Int64
            help = "Umi place after spliting by 'umi-string'. 1 based array."
            default = 2
            dest_name = "umi_p"

        "--verbose", "-v"
            help = "Print log messages to STDOUT."
            action = :store_true
            dest_name = "verb"

        "--use-grouped-reads", "-g"
            help = ("Force to use already grouped reads by UMI (using `umi_tools group`).
                    If so, provide which tag to use with --group-tag option.")
            action = :store_true
            dest_name = "use_grouped"

        "--group-tag", "-t"
            arg_type = String
            help = "Which tag should be used if --use-grouped-reads option is selected."
            default = "BX"
            dest_name = "group_tag"

        "--report-variants", "-R"
            help = "If to report sets of mutations."
            action = :store_true
            dest_name = "rep_vars"

        "--gene-conde", "-C"
            arg_type = Int64
            help = "Gene code from BioSequences.ncbi_trans_table"
            default = 11
            dest_name = "gene_code"

        "--translate", "-T"
            help = ("Set to translte genes. --genes-csv is needed.")
            action = :store_true
            dest_name = "translate"

        "--gene-csv", "-G"
            arg_type = String
            help = "Gene csv file if translate. Header: Chromosome,Start,End,Gene,Strand."
            dest_name = "gene_csv"

        "--algorithm", "-a"
            arg_type = String
            help = (
                "Algorithm for calculating fidelity. Available options are:
                umi - group reads to families/clusters by unique molecular identifiers.
                position - group reads to families/clusters by mapping positions.
                paired - do not group but determine confident error only if found in both paired end reads.
                    This means that insert of a library should be small to get as many overlaping/read-through as possible.
                naive - do not group and do not filter by anything.
                ")
            default = "umi"
            dest_name = "algo"

        "--frequency-in-cluster", "-F"
            arg_type = Float64
            help = "Mutation frequency in cluster. If larger, mutation passes."
            default = 1.0
            dest_name = "freq_in_clust"

        "--min-cluster-counts", "-L"
            arg_type = Int64
            help = "Minimum count of clusters for the mutation to be reported to a single_mutaitons file."
            default = 1
            dest_name = "min_cluster_counts"
        "--min-mutation-coverage", "-M"
            arg_type = Int64
            help = "Minimum confident coverage for the mutation to be reported to a single_mutaitons file."
            default = 1
            dest_name = "min_mutation_coverage"
    end

    parsed_args = parse_args(arg_parse_settings)

    file = parsed_args["bam"]
    alignQ = parsed_args["alignQ"]
    mut_n_read = parsed_args["mut_n_read"]
    q_limit = parsed_args["q_limit"]
    umi_s = parsed_args["umi_s"]
    umi_p = parsed_args["umi_p"]
    use_grouped = parsed_args["use_grouped"]
    group_tag = parsed_args["group_tag"]
    verbose = parsed_args["verb"]
    reference = parsed_args["ref"]
    clust_size = parsed_args["clust_size"]
    ref_error = parsed_args["ref_error"]
    output = parsed_args["out"]
    algo = parsed_args["algo"]
    freq_in_clust = parsed_args["freq_in_clust"]
    min_cluster_counts = parsed_args["min_cluster_counts"]
    min_mutation_coverage = parsed_args["min_mutation_coverage"]

    if !(algo in ["umi", "position", "paired", "naive"])
        println(STDERR, "ERROR: Selected --algorithm does not exist: $algo")
        println(STDERR, "ERROR: Options are: umi position paired naive")
        exit(1)
    end

    referenceSruct = ReferenceSequences(reference)
    cvCoverage = Coverage(referenceSruct)
    coverage_counter!(cvCoverage, file; mappingQuality=alignQ)

    clust_coverages = Coverage(referenceSruct)

    if algo == "paired" || algo == "naive"
        clust_size = 1
    end

    (vars,
    mut_coverages,
    raw_bases,
    total_reads,
    filtered_reads,
    filtered_R1,
    filtered_R2,
    error_distr_per_read
    ) = collect_variants(file, clust_coverages, alignQ, mut_n_read, q_limit,
                         umi_s, umi_p, use_grouped, group_tag; verbose=verbose, algo=algo)

    (passing_vars,
    err_distr_per_cluster,
    clust_size_distr,
    found_mutations,
    found_insertions,
    found_deletions,
    found_ns,
    total_clusters,
    conf_bases,
    clusts_with_vars,
    passing_mutations,
    total_bases,
    conf_coverage
    ) = filter_variants(vars, mut_coverages,
                           cvCoverage, referenceSruct,
                           clust_size, ref_error;
                           verbose=verbose, algo=algo,
                           freq_in_clust=freq_in_clust)

    df = DataFrame(
        Sample=String[],
        TotalErrors=Int64[],
        RawBases=Int64[],
        TotalClusterBases=Int64[],
        ConfidentBases=Int64[],
        Frequency=Float64[],
        MeanClustSize=Float64[],
        MinClusterSize=Int64[],
        MaxClusterSize=Int64[],
        SumClusterSize=Int64[],
        Mutations=Int64[],
        Insertions=Int64[],
        Deletions=Int64[],
        Ns=Int64[],
        TotalClusters=Int64[],
        ClustersWithVars=Int64[],
        TotalReads=Int64[],
        FilteredReads=Int64[],
        FilteredR1=Int64[],
        FilteredR2=Int64[]
    )

    sample_name = split(join(split(file,".")[1:end-1], "_"),"/")[end]
    push!(df, [
            sample_name,
            passing_vars,
            raw_bases,
            total_bases,
            conf_bases,
            passing_vars/conf_bases,
            mean(clust_size_distr),
            minimum(clust_size_distr),
            maximum(clust_size_distr),
            sum(clust_size_distr),
            found_mutations,
            found_insertions,
            found_deletions,
            found_ns,
            total_clusters,
            clusts_with_vars,
            total_reads,
            filtered_reads,
            filtered_R1,
            filtered_R2
        ]
    )

    CSV.write(output, df, delim=',')
    outfile = replace(output, ".csv", "")

    df = DataFrame(
        Sample=String[],
        Chromosome=String[],
        ClustOfVariants=String[],
        Count=Int64[],
        Coverage=Float64[],
        ConfidentCoverage=Float64[]
    )

    write_array(outfile*"_distributions.csv", err_distr_per_cluster)
    write_array(outfile*"_size_distribution.csv", clust_size_distr)
    write_array(outfile*"_error_distributions_per_read.csv", error_distr_per_read; new_line=false)

    transitions, transversions, types = transvers(passing_mutations)

    write_array(outfile*"_transitversions.csv", ["Sample","Transitions","Transversions"]; new_line=false)
    write_array(outfile*"_transitversions.csv", [sample_name, transitions, transversions]; new_line=false, meth="a")
    header = collect(keys(types))
    values = [string(types[i]) for i in header]
    write_array(outfile*"_transvers_types.csv", unshift!(header, "Sample"); new_line=false)
    write_array(outfile*"_transvers_types.csv", unshift!(values, sample_name); new_line=false, meth="a")

    if parsed_args["rep_vars"]
        mapped = mapp_mutations(passing_mutations)

        for mutations in keys(mapped)
            # average coverage
            cov = 0
            conf_cov = 0
            ct = 0
            clust_of_variants = ""
            chr = ""

            for mutation in mutations
                chr = mutation.chrom
                try
                    conf_cov += sum(conf_coverage.coverage[chr][mutation.mStart:mutation.mEnd])
                    cov += sum(clust_coverages.coverage[chr][mutation.mStart:mutation.mEnd])
                    ct += 1
                catch error
                    println(STDERR, "ERROR: $error")
                    continue
                end

                clust_of_variants = clust_of_variants*";"*stringify(mutation)
            end

            push!(df, [
                    sample_name,
                    chr,
                    clust_of_variants[2:end],
                    mapped[mutations],
                    round(cov/ct, 0),
                    round(conf_cov/ct, 0),
                ]
            )
        end

        CSV.write(outfile*"_passing_mutations.csv", df, delim=',')

        df = DataFrame(
            Sample=String[],
            Chromosome=String[],
            ClustOfVariants=String[],
            Count=Int64[],
            Coverage=Float64[],
            ConfidentCoverage=Float64[]
        )
        mapped = mapp_mutations(passing_mutations; cluster=false)

        for mutation in keys(mapped)
            # average coverage
            cov = 0
            conf_cov = 0
            chr = mutation.chrom

            try
                conf_cov += sum(conf_coverage.coverage[chr][mutation.mStart:mutation.mEnd])
                cov += sum(clust_coverages.coverage[chr][mutation.mStart:mutation.mEnd])
            catch error
                println(STDERR, "ERROR: $error")
                continue
            end

            if mapped[mutation] >= min_cluster_counts && conf_cov >= min_mutation_coverage
                push!(df, [
                        sample_name,
                        chr,
                        stringify(mutation),
                        mapped[mutation],
                        cov,
                        conf_cov,
                    ]
                )
            end
        end

        CSV.write(outfile*"_passing_single_mutations.csv", df, delim=',')

    end

    # Tranlate mutations

    if parsed_args["translate"]
        if verbose
            println(STDOUT, "INFO: translating passed mutations. Using ", parsed_args["gene_csv"])
        end

        if parsed_args["gene_csv"] == ""
            println(STDERR, "ERROR: gene bed file was not provided.")
            println(STDERR, "ERROR: Provide gene bed file if you want to translate mutations.")
            exit(1)
        end
        aa_substitutions, frame_shifts = translate_mutations(
            referenceSruct, parsed_args["gene_csv"],
            passing_mutations;
            genetic_code=parsed_args["gene_code"],
            verbose=verbose
            )

        df = DataFrame(
            Sample=String[],
            Chromosome=String[],
            ClustOfVariants=String[],
            Count=Int64[],
            Coverage=Float64[],
            ConfidentCoverage=Float64[]
        )

        for chr in keys(aa_substitutions)
            mapped = countmap(aa_substitutions[chr])
            for mutations in keys(mapped)
                # average coverage
                conf_cov = 0
                cov = 0
                ct = 0
                clust_of_variants = ""

                for mutation in mutations
                    position = mutation[2:end-1]
                    position = parse(Int64, String(position))
                    position = position * 3
                    try
                        conf_cov += sum(conf_coverage.coverage[chr][position-2:position])
                        cov += sum(clust_coverages.coverage[chr][position-2:position])
                        ct += 3
                    catch error
                        println(STDERR, "ERROR: $error")
                        continue
                    end

                    clust_of_variants = clust_of_variants*";"*mutation
                end

                if mapped[mutations] >= min_cluster_counts && round(conf_cov/ct, 0) >= min_mutation_coverage
                    push!(df, [
                            sample_name,
                            chr,
                            clust_of_variants[2:end],
                            mapped[mutations],
                            round(cov/ct, 0),
                            round(conf_cov/ct, 0)
                        ]
                    )
                end
            end
        end

        CSV.write(outfile*"_passing_aa_mutations.csv", df, delim=',')
        df = DataFrame(Sample=String[], FrameShifts=Int64[])
        push!(df, [sample_name frame_shifts])
        CSV.write(outfile*"_passing_frameshifts.csv", df, delim=',')

        df = DataFrame(
            Sample=String[],
            Chromosome=String[],
            Variants=String[],
            Count=Int64[],
            Coverage=Float64[],
            ConfidentCoverage=Float64[]
        )

        for chr in keys(aa_substitutions)
            mapped = mapp_single_aa_subs(aa_substitutions[chr])
            for mutation in keys(mapped)
                # average coverage
                conf_cov = 0
                cov = 0
                ct = 0
                clust_of_variants = ""

                position = mutation[2:end-1]
                position = parse(Int64, String(position))
                position = position * 3
                try
                    conf_cov += sum(conf_coverage.coverage[chr][position-2:position])
                    cov += sum(clust_coverages.coverage[chr][position-2:position])
                    ct += 3
                catch error
                    println(STDERR, "ERROR: $error")
                    continue
                end

                if mapped[mutation] >= min_cluster_counts && round(conf_cov/ct, 0) >= min_mutation_coverage
                    push!(df, [
                            sample_name,
                            chr,
                            mutation,
                            mapped[mutation],
                            round(cov/ct, 0),
                            round(conf_cov/ct, 0)
                        ]
                    )
                end
            end
        end
        CSV.write(outfile*"_passing_single_aa_mutations.csv", df, delim=',')
    end
end


function write_array(flname::String, arr; new_line::Bool=true, meth::String="w")

    out = open(flname, meth)

    if new_line
        for i in arr
            write(out, string(i)*"\n")
        end
    else
        ct = 0
        for i in arr
            ct += 1
            if ct == 1
                write(out, string(i))
            else
                write(out, ","*string(i))
            end
        end
        write(out, "\n")
    end

    close(out)
end

main(ARGS)
