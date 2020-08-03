#!/usr/bin/env julia

# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

__precompile__()

module FidelityUMI

export
    ReferenceSequences,
    Coverage,
    Mutation,
    filter_variants,
    call_variants,
    check_operations,
    check_quality,
    collect_variants,
    count_mut_coverages!,
    coverage_counter!,
    dedup_mutations,
    get_read_shift,
    populate_variants!,
    translate_mutations,
    compare_seq,
    mapp_mutations,
    parse_chr,
    stringify,
    conf_coverage!,
    transvers,
    MutationsCluster,
    mapp_single_aa_subs


import BioAlignments: BioAlignments, BAM, Operation
import BioAlignments: OP_MATCH, OP_INSERT, OP_DELETE, OP_SKIP
import BioAlignments: OP_SOFT_CLIP, OP_HARD_CLIP, OP_PAD
import BioAlignments: OP_SEQ_MATCH, OP_SEQ_MISMATCH, OP_BACK, OP_START
import BioSequences: DNASequence, FASTA, ReferenceSequence, BioSequences
import BioSequences: sequence, seqname, RNASequence, translate
import AutoHashEquals: @auto_hash_equals
import CSV: CSV, eachrow
import DataFrames: DataFrame


type  ReferenceSequences
    file::String
    chromosomes::Array{String,1}
    chrom_lengths::Dict{String,Int64}
    sequences::Dict{String,BioSequences.ReferenceSequence}

    function ReferenceSequences(file::String)

        sequences::Dict{String,BioSequences.ReferenceSequence} = Dict()
        chromosomes::Array{String,1} = []
        chrom_lengths::Dict{String,Int64} = Dict()
        ext = split(file, ".")[end]

        if ext == "fasta" || ext == "fa" || ext == "fna"
            reader = FASTA.Reader(open(file,"r"))

        elseif ext == "2bit"
            reader = TwoBit.Reader(open(file, "r"))

        else
            println(STDERR, "ERROR: Reference FILE extention is not recognised: $ext")
            ext(1)
        end

        for record in reader
            sequences[seqname(record)] = convert(DNASequence,replace(string(sequence(record)),r"[^ATGCN]","N"))
            push!(chromosomes,seqname(record))
            chrom_lengths[seqname(record)] = length(sequences[seqname(record)])
        end

        return new(file, chromosomes, chrom_lengths, sequences)::ReferenceSequences
    end
end


type  Coverage
    chromosomes::Array{String,1}
    chrom_lengths::Dict{String,Int64}
    coverage::Dict{String,Array{UInt32,1}}

    function Coverage(reference::ReferenceSequences)
        chromosomes::Array{String,1} = reference.chromosomes
        chrom_lengths::Dict{String,Int64} = reference.chrom_lengths
        coverage::Dict{String,Array{UInt32,1}} = Dict()

        for chr in chromosomes
            coverage[chr] = zeros(UInt32, chrom_lengths[chr])
        end

        return new(chromosomes, chrom_lengths, coverage)
    end
end


@auto_hash_equals immutable Mutation
    class::Char
    chrom::String
    mStart::Int64
    mEnd::Int64
    reference::String
    variant::String
end


@auto_hash_equals mutable struct MutationsCluster
    mutation_count::Dict{Mutation,Int64}
    cluster_size::Int64
    R1_length::Int64
    R2_length::Int64
    total_mutations::Int64
    R1_coverage::Array{UInt16,1}
    R2_coverage::Array{UInt16,1}
    c_start::Int64
    c_end::Int64
end


function coverage_counter!(Coverage::Coverage, bamF::String; mappingQuality::Int64=0)
    reader = open(BAM.Reader, bamF)
    record = BAM.Record()
    ct::Int64 = 0
    ct2::Int64 = 0
    while !eof(reader)
        ct += 1
        read!(reader, record)

        if BAM.ismapped(record)
            refName::String = BAM.refname(record)
            flag = BAM.flag(record)
            mpq = BAM.mappingquality(record)
            # these conditions should match collect_variants() conditions
            if mpq >= mappingQuality && (flag&4 == 0) && (flag&256 == 0) && (flag&2048 == 0) && (refName in Coverage.chromosomes)   # filter out only mapped and primary alignments
                Coverage.coverage[refName][BAM.position(record):BAM.rightposition(record)]+=UInt32(1)
                ct2 += 1
            end
        end
    end
    close(reader)
end


function call_variants(record::BioAlignments.BAM.Record; quality::Int64=0, verbose::Bool=false)

    mutations::Array{Mutation,1} = []
    leftposition = BAM.position(record)
    md_tag = uppercase(record["MD"]::String)
    record_cigar = BAM.cigar(record)
    rfname = BAM.refname(record)
    read_seq = BAM.sequence(record)
    a_operations, o_lengths = BAM.cigar_rle(record)
    variants, md_positions = split_md_tag(md_tag)
    qualities = BAM.quality(record)

    if !(check_operations(a_operations))
        println(STDOUT, "INFO: CALLING VARIANTS, IN RECORD:")
        println(STDOUT, record)
        println()
        return Mutation[]
    end

    if verbose
        println(STDOUT, "INFO: CALLING VARIANTS, MD_TAG: ", md_tag)
        println(STDOUT, "INFO: CALLING VARIANTS, PARSED VARIANTS AND POSITIONS: ", variants, ":", md_positions)
        println(STDOUT, "INFO: CALLING VARIANTS, CIGAR RECORD: ", record_cigar)
        println(STDOUT, "INFO: CALLING VARIANTS, ALIGNMENT OPERATIONS: ", a_operations)
        println(STDOUT, "INFO: CALLING VARIANTS, LENGTHS OF OPERATIONS: ", o_lengths)
    end

    # Find variants and deletions in the read.
    ct = 0
    for i in variants
        ct += 1
        ref_shift = get_ref_shift(variants, md_positions, ct)
        read_shift = get_read_shift(a_operations, o_lengths, variants, md_positions, ct)

        if verbose
            println(STDOUT, "INFO: CALL VARIANTS, REF_SHIFT: ", ref_shift)
            println(STDOUT, "INFO: CALL VARIANTS, READ_SHIFT: ", read_shift)
        end

        if i[1] == '^'
                class = 'D'
            ix = 2
            s2 = length(i[2:end]) - 1
        else
            class = 'M'
            ix = 1
            s2 = 0
        end

        var_start_pos = leftposition + ref_shift - s2 - 1
        var_end_pos = leftposition + ref_shift - 1

        if verbose
            println(STDOUT, "INFO: CALL VARIANTS, READ LEFT MAPPING POSITION: ", leftposition)
            println(STDOUT, "INFO: CALL VARIANTS, CALCULATED SHIFT in REFERENCE: ", ref_shift)
            println(STDOUT, "INFO: CALL VARIANTS, VARIANT START POSITION IN REFERENCE: ", var_start_pos)
            println(STDOUT, "INFO: CALL VARIANTS, VARIANT END POSITION IN REFERENCE: ", var_end_pos)
            println(STDOUT, "INFO: CALL VARIANTS, QUALITY MAPPING START POSITION: ", read_shift)
            println(STDOUT, "INFO: CALL VARIANTS, QUALITY MAPPING end POSITION: ", read_shift + s2)
        end

        if i[1] == '^'
            if qualities[read_shift] >= quality
                push!(mutations, Mutation(class, rfname, var_start_pos, var_end_pos, i[ix:end], "*"))
            else
                if verbose
                    println(STDOUT, "SKIPED BY QUALITY: ", mean(qualities[read_shift:read_shift + s2]))
                end
            end
        else
            if mean(qualities[read_shift:read_shift]) >= quality
                v = read_seq[read_shift:read_shift]
                push!(mutations, Mutation(class, rfname, var_start_pos, var_end_pos, i[ix:end], convert(String, v)))
            else
                if verbose
                    println(STDOUT, "INFO: SKIPED BY QUALITY: ", mean(qualities[read_shift:read_shift]), " ", qualities[read_shift:read_shift])
                end
            end
        end
    end


    # find insertions in the read
    ct = 0
    read_pos = 0
    ref_pos = 0

    for i in a_operations
        ct += 1

        if i == BioAlignments.Operation(OP_MATCH)
            ref_pos += o_lengths[ct]
            read_pos += o_lengths[ct]

        elseif i == BioAlignments.Operation(OP_DELETE)
            ref_pos += o_lengths[ct]

        elseif i == BioAlignments.Operation(OP_SOFT_CLIP)
            read_pos += o_lengths[ct]

        elseif i == BioAlignments.Operation(OP_INSERT)
            in_start = read_pos + 1
            read_pos += o_lengths[ct]
            in_end = read_pos
            v = read_seq[in_start:in_end]
            var_start_pos = ref_pos + leftposition - 1
            var_end_pos = var_start_pos + 1

            if mean(qualities[in_start:in_end]) >= quality
                push!(mutations, Mutation('I', rfname, var_start_pos, var_end_pos, "-", v))
            end
        end
    end

    if verbose
        println(STDOUT, "INFO: VARIANTS CALCULATED: ", mutations)
        println(STDOUT, "INFO: FROM THIS RECORD: ", record)
        println()
    end

    return mutations
end


function check_operations(op::Array{BioAlignments.Operation,1})

    supported_operations = [BioAlignments.Operation(OP_MATCH),
                            BioAlignments.Operation(OP_INSERT),
                            BioAlignments.Operation(OP_SOFT_CLIP),
                            BioAlignments.Operation(OP_DELETE),
                            BioAlignments.Operation(OP_HARD_CLIP)]

    for i in op
        if !(i in supported_operations)
            println(STDOUT, "WARNING: FOUND NOT SUPPORTED ALIGMENT OPERATION: $i")
            return false
        end
    end

    return true
end


function split_md_tag(md::String)
    md_split::Array{String, 1} = split(replace(md, r"[\\^]*[ACGT]+", x -> "_"*x*"_"), "_")
    variants = String[]
    positions = Int64[]

    for i in md_split
        # finding numbers
        if ismatch(r"[0-9]+", i)
             pos = parse(Int64, i)
             push!(positions, pos)
        else
            push!(variants, i)
        end
    end

    return variants, positions
end


function get_read_shift(cigar_o::Array{BioAlignments.Operation,1},
                        lengths::Array{Int64,1},
                        variants::Array{String,1},
                        md_positions::Array{Int64,1},
                        pos::Int64
                        )
    shift = 0
    op_shitfs = 0
    ct = 0

    for i in variants
        ct += 1
        if pos >= ct && i[1] != '^'
            shift += md_positions[ct] + 1
            op_shitfs += md_positions[ct] + 1
        elseif pos >= ct
            shift += md_positions[ct]
            op_shitfs += md_positions[ct]
        end
    end

    ct = 0
    op_matches = 0
    for i in cigar_o
        ct += 1
        if i == BioAlignments.Operation(OP_MATCH)
            op_matches += lengths[ct]

        elseif i == BioAlignments.Operation(OP_SOFT_CLIP)
            shift += lengths[ct]

        elseif i == BioAlignments.Operation(OP_INSERT)
            shift += lengths[ct]
        end

        if  op_shitfs < op_matches
            break
        end
    end

    return shift
end


function get_ref_shift(variants::Array{String,1}, md_positions::Array{Int64,1}, pos::Int64)

    ct = 0
    shift = 0
    for i in md_positions
        ct += 1
        if pos >= ct
            shift += i
            if variants[ct][1] == '^'
                shift += length(variants[ct][2:end])
            else
                shift += 1
            end
        end
    end

    return shift
end


function collect_variants(bamfile::String, Coverage::Coverage, mappingq::Int64, mut_n_read::Int64, q_limit::Int64,
                          umi_s::String, umi_p::Int64, use_grouped::Bool, group_tag::String;
                          verbose::Bool=false, algo::String="umi")

    if !(algo in ["umi", "position", "paired", "naive"])
        println(STDERR, "ERROR: Selected --algorithm does not exist: $algo")
        println(STDERR, "ERROR: Options are: umi position paired naive")
        exit(1)
    end

    # collect coverages for each mutation from all clusters
    # TODO: Change mutcoverage structure to remove chromosome as key since Mutation has chromosome and is unique.
    mut_coverages = Dict{String,Dict{Mutation,Int64}}()
    # key is unique ID based either on mapping position or UMI + mapping position
    vars = Dict{String,MutationsCluster}()

    error_distr = Int64[]
    read_distr = Int64[]
    reader = open(BAM.Reader, bamfile)
    record = BAM.Record()
    record2 = BAM.Record()
    total_bases = 0
    total_reads = 0
    filtered_reads = 0
    filtered_R1 = 0
    filtered_R2 = 0

    while !eof(reader)
        read!(reader, record)

        try
            read!(reader, record2)
        catch
            break
        end

        if BAM.ismapped(record) && BAM.ismapped(record2)
            passed1 = false
            passed2 = false
            flag1 = BAM.flag(record)
            flag2 = BAM.flag(record2)
            tempname::String = BAM.tempname(record)
            tempname2::String = BAM.tempname(record2)
            rfname::String = BAM.refname(record)
            rfname2::String = BAM.refname(record2)
            mpq = BAM.mappingquality(record)
            mpq2 = BAM.mappingquality(record2)
            rpos = BAM.rightposition(record)
            lpos = BAM.position(record)
            rpos2 = BAM.rightposition(record2)
            lpos2 = BAM.position(record2)
            seq_len = BAM.seqlength(record)
            seq_len2 = BAM.seqlength(record2)
            mutations::Array{Mutation,1} = call_variants(record; quality=q_limit, verbose=verbose)
            mutations2::Array{Mutation,1} = call_variants(record2; quality=q_limit, verbose=verbose)
            l1 = length(mutations) > 0
            l2 = length(mutations2) > 0

            # This section is for counting mutations at positions
            # Conditions must match coverage_counter!() conditions
            if  mpq >= mappingq && flag1&4 == 0 && flag1&256 == 0 && flag1&2048 == 0 && l1 <= mut_n_read
                total_reads += 1
                passed1 = true
                if l1
                    # count before deduplication of overlaped paired regions
                    count_mut_coverages!(mut_coverages, mutations, rfname; verbose=verbose)
                    error_dist_over_length!(error_distr, mutations, lpos, seq_len)
                    read_dist!(read_distr, seq_len)
                end
            end

            if mpq2 >= mappingq && flag2&4 == 0 && flag2&256 == 0 && flag2&2048 == 0 && l2 <= mut_n_read
                total_reads += 1
                passed2 = true
                if l2
                    # count before deduplication of overlaped paired regions
                    count_mut_coverages!(mut_coverages, mutations2, rfname2; verbose=verbose)
                    error_dist_over_length!(error_distr, mutations2, lpos2, seq_len2)
                    read_dist!(read_distr, seq_len2)
                end
            end

            # only if first is R1 and second is R2
            # only if both mapped and primary alignments
            if tempname == tempname2 && passed1 && passed2 && rfname == rfname2
                if use_grouped && algo == "umi"
                    try
                        if flag1&64 == 64
                            UMI = uppercase(record[group_tag]::String)
                        else
                            UMI = uppercase(record2[group_tag]::String)
                        end
                    catch e
                        println(STDERR, "WARNING: GROUP TAG WAS NOT FOUND: $group_tag")
                        println(STDERR, "TRYING UMI FROM READ NAME ....")
                        UMI = split(tempname, umi_s)[umi_p]
                    end
                elseif algo == "umi"
                    UMI = split(tempname, umi_s)[umi_p]
                end

                spos = minimum([lpos, lpos2, rpos, rpos2])
                epos = maximum([lpos, lpos2, rpos, rpos2])

                if algo == "umi"
                    key = rfname*":"*string(spos)*":"*string(epos)*":"*UMI
                elseif algo == "position"
                    key = rfname*":"*string(spos)*":"*string(epos)
                else
                    key = rfname*":"*string(spos)*":"*string(epos)*":"*tempname
                end

                if algo == "paired"
                    # duplicated mutations found in R1 and R2 overlaping reads
                    variants = intersect(mutations, mutations2)
                else
                    # unique mutations from R1 and R2
                    variants = union(mutations, mutations2)
                end

                # if array is not empty after base quality filtering.

                if flag1&16 == 0
                    seq_len = rpos - lpos + 1
                    seq_len2 = rpos2 - lpos2 + 1
                else
                    seq_len2 = rpos - lpos + 1
                    seq_len = rpos2 - lpos2 + 1
                end

                total_bases += seq_len
                total_bases += seq_len2

                if haskey(vars, key)
                    populate_variants!(vars, variants, key, seq_len,
                                       seq_len2, spos, epos;
                                       verbose=verbose)
                else
                    Coverage.coverage[rfname][lpos:rpos] += 1
                    Coverage.coverage[rfname][lpos2:rpos2] += 1
                    # length = epos - spos + 1
                    populate_variants!(vars, variants, key, seq_len,
                                       seq_len2, spos, epos;
                                       verbose=verbose)
                end
            else
                if tempname != tempname2
                    println(STDOUT, "WARNING: READ NAMES DOES NOT MATCH. Is bam sorted by name?")
                    filtered_reads += 2
                end

                if rfname2 != rfname
                    println(STDOUT, "WARNING: READS ARE ON DIFF CHROMOSOMES [$tempname]")
                    filtered_reads += 2
                end

                if !(passed1)
                    filtered_R1 += 1
                end

                if !(passed2)
                    filtered_R2 += 1
                end
            end
        end
    end

    if verbose
        println(STDOUT, "INFO: COVERAGES OF MUTATIONS: ")
        for i in mut_coverages
            println(STDOUT, i)
        end
        println(STDOUT, "INFO: END OF COVERAGES OF MUTATIONS: ")
    end

    error_distr = error_distr./read_distr

    return vars, mut_coverages, total_bases, total_reads, filtered_reads, filtered_R1, filtered_R2, error_distr
end


function count_mut_coverages!(mut_coverages::Dict{String,Dict{Mutation,Int64}},
                              mutations::Array{Mutation,1}, refname::String;
                              verbose::Bool=false
                              )

    if verbose
        println(STDOUT, "INFO: COUNT MUTATION COVERAGES, input: $mut_coverages")
    end


    for mutation in mutations
        if verbose
            println(STDOUT, "INFO: COUNT MUTATION COVERAGES, MUTATION: $mutation")
        end

        # collecting coverages of mutations per reference
        if haskey(mut_coverages, refname)
            if haskey(mut_coverages[refname], mutation)
                mut_coverages[refname][mutation] += 1
            else
                mut_coverages[refname][mutation] = 1
            end
        else
            mut_coverages[refname] = Dict(mutation => 1)
        end

        if verbose
            println(STDOUT, "INFO: COUNT MUTATION COVERAGES, MUTATION COVERAGE: ", mut_coverages[refname][mutation])
        end

    end

    return mut_coverages
end


function populate_variants!(vars::Dict{String,MutationsCluster},
                            mutations::Array{Mutation, 1},
                            id::String, seq_len::Int64,
                            seq_len2::Int64,
                            spos::Int64, epos::Int64;
                            verbose::Bool=false)

    if !(haskey(vars, id))
        if verbose
            println(STDOUT, "INFO: POPULATE VARIANTS, init new cluster w/o mutations: R1_len:$seq_len R2_len:$seq_len2 spos:$spos epos:$epos, id:$id")
        end
        vars[id] = MutationsCluster(
            Dict(),
            0,
            seq_len,
            seq_len2,
            0,
            zeros(UInt16, seq_len),
            zeros(UInt16, seq_len2),
            spos,
            epos
        )
    end

    for mutation in mutations
        if verbose
            println(STDOUT, "INFO: POPULATE VARIANTS, adding $mutation $seq_len $seq_len2 $spos $epos to cluster $id", vars[id] )
        end
        if haskey(vars[id].mutation_count, mutation)
            vars[id].mutation_count[mutation] += 1
        else
            vars[id].mutation_count[mutation] = 1
            vars[id].total_mutations += 1
        end
    end

    if vars[id].R1_length < seq_len
        if verbose
            println(STDOUT, "INFO: POPULATE VARIANTS, R1_len increased: ", vars[id].R1_length, " $seq_len")
            println(STDOUT, "INFO: POPULATE VARIANTS, extending R1_coverage: ", length(vars[id].R1_coverage))
        end
        append!(vars[id].R1_coverage, zeros(seq_len-vars[id].R1_length))
        if verbose
            println(STDOUT, "INFO: POPULATE VARIANTS, extended R1_coverage: ", length(vars[id].R1_coverage))
        end
        vars[id].R1_length = seq_len
    end

    if vars[id].R2_length < seq_len2
        if verbose
            println(STDOUT, "INFO: POPULATE VARIANTS, R2_len increased: ", vars[id].R1_length, " $seq_len")
            println(STDOUT, "INFO: POPULATE VARIANTS, extending R1_coverage: ", length(vars[id].R1_coverage))
        end
        append!(vars[id].R2_coverage, zeros(seq_len2-vars[id].R2_length))
        if verbose
            println(STDOUT, "INFO: POPULATE VARIANTS, extended R1_coverage: ", length(vars[id].R1_coverage))
        end
        vars[id].R2_length = seq_len2
    end

    if verbose
        println(STDOUT, "INFO: POPULATE VARIANTS, incrementing R1_coverage: ", vars[id].R1_coverage)
    end

    vars[id].R1_coverage[1:seq_len] += 1

    if verbose
        println(STDOUT, "INFO: POPULATE VARIANTS, incremented R1_coverage: ", vars[id].R1_coverage)
        println(STDOUT, "INFO: POPULATE VARIANTS, incrementing R2_coverage: ", vars[id].R2_coverage)
    end

    vars[id].R2_coverage[1:seq_len2] += 1

    if verbose
        println(STDOUT, "INFO: POPULATE VARIANTS, incremented R2_coverage: ", vars[id].R2_coverage)
        println(STDOUT, "INFO: POPULATE VARIANTS, incrementing cluster size: ", vars[id].cluster_size)
    end

    vars[id].cluster_size += 1

    if verbose
        println(STDOUT, "INFO: POPULATE VARIANTS, incremented cluster size: ", vars[id].cluster_size)
    end

    return vars
end

function filter_variants(vars::Dict{String,MutationsCluster},
                     mut_coverages::Dict{String,Dict{Mutation,Int64}},
                     cvCoverage::Coverage, referenceSruct::ReferenceSequences,
                     cluster_size::Int64, ref_error::Float64; algo::String="umi",
                     freq_in_clust::Float64=1.0, verbose::Bool=false)

    if !(algo in ["umi", "position", "paired", "naive"])
        println(STDERR, "ERROR: Selected --algorithm does not exist: $algo")
        println(STDERR, "ERROR: Options are: umi position paired naive")
        exit(1)
    end

    true_var_count = 0
    err_distr = Float32[]
    clust_size_distr = Int64[]
    clusts_with_vars = 0
    total_clusters = 0
    found_insertions = 0
    found_deletions = 0
    found_mutations = 0
    found_ns = 0
    passing_vars = Dict{String,MutationsCluster}()
    confident_bp = 0
    total_bp = 0
    conf_coverage = Coverage(referenceSruct)

    for id in keys(vars)
        # parse key and get position of mutation
        chromosome = parse_chr(id)
        read_position = parse(Int64, split(id, ":")[2])
        found_var = false
        total_clusters += 1

        if verbose
            println("INFO: CALCULATE VARIANTS, CLUSTER ID: ", id)
            println("INFO: CALCULATE VARIANTS, CLUSTER: ", vars[id])
        end


        push!(clust_size_distr, vars[id].cluster_size)
        if vars[id].cluster_size >= cluster_size

            # ad cluster to confident coverage
            start = vars[id].c_start
            end_ = vars[id].c_start + vars[id].R1_length - 1
            try
                conf_coverage.coverage[chromosome][start:end_] += 1
            catch error
                println(STDERR, "ERROR: $error")
                println(STDERR, "ERROR: ", vars[id])
            end
            end_ = vars[id].c_end
            start = vars[id].c_end - vars[id].R2_length + 1
            try
                conf_coverage.coverage[chromosome][start:end_] += 1
            catch error
                println(STDERR, "ERROR: $error")
                println(STDERR, "ERROR: ", vars[id])
            end

            for variant in keys(vars[id].mutation_count)
                var_count = vars[id].mutation_count[variant]
                spos = variant.mStart
                epos = variant.mEnd

                # Filtering non confident variants in the cluster
                if var_count/vars[id].cluster_size >= freq_in_clust
                    mut_cov = mut_coverages[chromosome][variant]
                    coverage = mean(cvCoverage.coverage[chromosome][spos:epos])

                    # checking if mutation coverage is not equal to coverage
                    if mut_cov/coverage < ref_error

                        if verbose
                            println("INFO: CALCULATE VARIANTS, CLUSTER ID: ", id)
                            println("INFO: CALCULATE VARIANTS, READ_POSITION: ", read_position)
                            println("INFO: CALCULATE VARIANTS, COVERAGE: ", coverage)
                            println("INFO: CALCULATE VARIANTS, MUT_COVERAGE: ", mut_cov)
                            println("INFO: CALCULATE VARIANTS, ERROR_COUNT: ", var_count)
                        end

                        found_var = true
                        true_var_count += 1

                        if variant.class == 'M'
                            found_mutations += 1
                        elseif variant.class == 'I'
                            found_insertions += 1
                        elseif variant.class == 'D'
                            found_deletions += 1
                        elseif variant.class == 'N'
                            found_ns += 1
                        end

                        if haskey(passing_vars, id)
                            passing_vars[id].mutation_count[variant] = 1
                            passing_vars[id].total_mutations += 1
                        else
                            passing_vars[id] = MutationsCluster(
                            Dict(variant => 1),
                            vars[id].cluster_size,
                            vars[id].R1_length,
                            vars[id].R2_length,
                            1,
                            vars[id].R1_coverage,
                            vars[id].R2_coverage,
                            vars[id].c_start,
                            vars[id].c_end
                            )
                        end

                        if verbose
                            println(STDOUT, "INFO: CALCULATE VARIANTS, PASSING MUTATION: $variant")
                        end
                    end

                else
                    if verbose
                        println("INFO: CALCULATE VARIANTS, CLUSTER ID: ", id)
                        println("INFO: CALCULATE VARIANTS, SKIPPING VAR: ", variant)
                    end

                end

                # collect error frequencies for distribution ploting.
                if var_count != 0
                    push!(err_distr, round(var_count/vars[id].cluster_size, 2))
                end
            end

            if algo == "umi" || algo == "position" || algo == "naive"
                # adding only confident bases, equal to cluster coverage.
                conf_bp = count(x -> x >= vars[id].cluster_size, vars[id].R1_coverage)
                t_bp = conf_bp * vars[id].cluster_size
            elseif algo == "paired"
                # adding only confident bases from overlaping pairs
                if vars[id].c_start + vars[id].R1_length > vars[id].c_end - vars[id].R2_length && vars[id].total_mutations > 0
                    conf_bp = vars[id].c_start + vars[id].R1_length - (vars[id].c_end - vars[id].R2_length) + 1
                    t_bp = vars[id].R1_length + vars[id].R2_length
                else
                    conf_bp = 0
                    t_bp = 0
                end
            end

            if verbose
                println(STDOUT, "INFO: CALCULATE VARIANTS, counted confident bp for R1 $id: $conf_bp")
                println(STDOUT, "INFO: CALCULATE VARIANTS, counted total cluster bp for R1 $id: $t_bp")
            end

            confident_bp += conf_bp
            total_bp += t_bp

            if algo == "umi" || algo == "position" || algo == "naive"
                conf_bp = count(x -> x >= vars[id].cluster_size, vars[id].R2_coverage)
                t_bp = conf_bp * vars[id].cluster_size

                if verbose
                    println(STDOUT, "INFO: CALCULATE VARIANTS, counted confident bp for R2 $id: $conf_bp")
                    println(STDOUT, "INFO: CALCULATE VARIANTS, counted total cluster bp for R2 $id: $t_bp")
                end

                confident_bp += conf_bp
                total_bp += t_bp
            end

        else
            if verbose
                println("INFO: CALCULATE VARIANTS, CLUSTER ID: ", id)
                println("INFO: CALCULATE VARIANTS, SKIPPING CLUSTER: ", vars[id])
            end
        end

        if found_var
            clusts_with_vars += 1
        end
    end

    return (
        true_var_count,
        err_distr,
        clust_size_distr,
        found_mutations,
        found_insertions,
        found_deletions,
        found_ns,
        total_clusters,
        confident_bp,
        clusts_with_vars,
        passing_vars,
        total_bp,
        conf_coverage
    )
end


function error_dist_over_length!(error_distr::Array{Int64, 1}, mutations::Array{Mutation,1}, lpos::Int64, read_length::Int64)

    if read_length > length(error_distr)
        for i in 1:(read_length-length(error_distr))
            push!(error_distr, 0)
        end
    end

    # find if any deletions

    l = 0
    # TODO: mutations should be sorted.
    for mutation in mutations
        if mutation.class == 'D'
            l += length(mutation.reference)
        end
        position_in_read = mutation.mStart - lpos - l + 1
        error_distr[position_in_read] += 1
    end

    return error_distr
end


function read_dist!(read_distr::Array{Int64, 1}, read_length::Int64)

    if read_length > length(read_distr)
        for i in 1:(read_length-length(read_distr))
            push!(read_distr, 0)
        end
    end

    read_distr .+= 1

    return read_distr
end


"""
    Translates collection of mutations to collection of aminoacid subtitutions.
    For genetic codes check BioSequences.ncbi_trans_table.
"""
function translate_mutations(referenceSruct::ReferenceSequences, genesCSV::String,
                             vars::Dict{String,MutationsCluster};
                             genetic_code::Int64=11, verbose::Bool=false)

    # get AA sequence from reference
    # loop mutations, get codon, translate.

    if verbose
        println(STDOUT, "INFO: TRANSLATE MUTATIONS, Loading BED file. $genesCSV")
    end
    df = CSV.read(genesCSV, delim=',')
    proteins = BioSequences.AminoAcidSequence[]
    genes_seqs = BioSequences.BioSequence{BioSequences.DNAAlphabet{4}}[]
    chromosomes = String[]
    gene_sidx = Int64[]
    gene_eidx = Int64[]
    strands = String[]
    names = String[]
    frame_shifts = 0
    # chr => variant_clusters[variant_cluster[variants]]
    aa_substitutions = Dict{String,Array{Array{String,1},1}}()

    if verbose
        println(STDOUT, "INFO: TRANSLATE MUTATIONS, Parsing records.")
    end

    for row in eachrow(df)
        if verbose
            println(STDOUT, "INFO: TRANSLATE MUTATIONS, Parsing recod: $row")
        end

        chr = row[1]
        strand = row[5]
        gene_start = row[2]
        gene_end = row[3]
        gene_name = row[4]

        if verbose
            println(STDOUT, "INFO: TRANSLATE MUTATIONS, strand: $strand")
            println(STDOUT, "INFO: TRANSLATE MUTATIONS, chromosome: $chr")
            println(STDOUT, "INFO: TRANSLATE MUTATIONS, gene start: $gene_start")
            println(STDOUT, "INFO: TRANSLATE MUTATIONS, gene end: $gene_end")
            println(STDOUT, "INFO: TRANSLATE MUTATIONS, gene name: $gene_name")
            println(STDOUT, "INFO: TRANSLATE MUTATIONS, loading gene sequence from reference.")
        end

        dna_sequence = DNASequence(referenceSruct.sequences[chr][gene_start:gene_end])
        if verbose
            println(STDOUT, dna_sequence)
            println(STDOUT, "INFO: TRANSLATE MUTATIONS, loading gene sequence from reference. [DONE]")
            println(STDOUT, "INFO: TRANSLATE MUTATIONS, convert DNA sequence to RNA.")
        end

        push!(genes_seqs, dna_sequence)

        if strand == "-"
            if verbose
                println(STDOUT, "INFO: TRANSLATE MUTATIONS, geting reverse complement.")
            end
            BioSequences.reverse_complement!(dna_sequence)
            if verbose
                println(STDOUT, "INFO: TRANSLATE MUTATIONS, geting reverse complement. [DONE]")
            end
        end

        rna_sequence = convert(RNASequence, dna_sequence)
        if verbose
            println(STDOUT, rna_sequence)
            println(STDOUT, "INFO: TRANSLATE MUTATIONS, convert DNA sequence to RNA. [DONE]")
        end

        push!(gene_sidx, gene_start)
        push!(gene_eidx, gene_end)
        push!(chromosomes, chr)
        push!(strands, strand)
        push!(names, gene_name)

        if verbose
            println(STDOUT, "INFO: TRANSLATE MUTATIONS, translating to protein sequence.")
        end
        protein = translate(rna_sequence, code=BioSequences.ncbi_trans_table[genetic_code], allow_ambiguous_codons=true)
        if verbose
            println(STDOUT, "INFO: TRANSLATE MUTATIONS, CHR: $chr; STRAND: $strand; START: $gene_start; END: $gene_end; PROTEIN: $protein")
        end


        push!(proteins, protein)
    end

    if verbose
        println(STDOUT, "INFO: TRANSLATE MUTATIONS, slelected genes chromosomes: $chromosomes")
        println(STDOUT, "INFO: TRANSLATE MUTATIONS, slelected genes starts: $gene_sidx")
        println(STDOUT, "INFO: TRANSLATE MUTATIONS, slelected genes ends: $gene_eidx")
        println(STDOUT, "INFO: TRANSLATE MUTATIONS, slelected genes strands: $strands")
    end

    # mutaion is Chromosome:ReferencePositionMutation
    for id in keys(vars)
        chr = parse_chr(id)

        if verbose
            println(STDOUT, "INFO: TRANSLATE MUTATIONS, cluster of mutations at $chr: ", vars[id])
        end
        perm_seq = DNASequence()
        gene_idx = 0
        found_first = false
        found_gene = false
        found_indel = false
        indel_ct = 0

        for var in keys(vars[id].mutation_count)
            mut_pos = var.mStart
            # search for gene matching cluster range
            if verbose
                println(STDOUT, "INFO: TRANSLATE MUTATIONS, Searching if var $var is matching any gene from BED.")
            end

            chr_idxs = find(x -> x == chr, chromosomes)
            gene_starts = find(x -> x <= mut_pos, gene_sidx)
            gene_ends = find(x -> x >= mut_pos, gene_eidx)
            common1 = intersect(chr_idxs, gene_starts)
            gene_idx2 = intersect(common1, gene_ends)

            if verbose
                println(STDOUT, "INFO: TRANSLATE MUTATIONS, Matching chromosomes indexes: $chr_idxs")
                println(STDOUT, "INFO: TRANSLATE MUTATIONS, Matching gene starts: $gene_starts")
                println(STDOUT, "INFO: TRANSLATE MUTATIONS, Matching gene ends indexes: $gene_ends")
                println(STDOUT, "INFO: TRANSLATE MUTATIONS, Matching both chromosomes and gene starts: $common1")
                println(STDOUT, "INFO: TRANSLATE MUTATIONS, Matching all: $gene_idx2")
            end

            if length(gene_idx2) < 1
                if verbose
                    println(STDOUT, "INFO: TRANSLATE MUTATIONS, $var is not in your gene list therefore will not be translated.")
                end
                continue
            elseif length(gene_idx2) > 1
                if verbose
                    println(STDOUT, "WARNING: TRANSLATE MUTATIONS, $var is matching more than one gene. Does ranges overlap? This will not be handled!")
                end
                continue
            else
                found_gene = true
                gene_idx = gene_idx2
                if verbose
                    println(STDOUT, "INFO: TRANSLATE MUTATIONS, $var found in ", names[gene_idx][1])
                end
            end

            # get gene sequence
            if !(found_first)
                perm_seq = copy(genes_seqs[gene_idx[1]])
                found_first = true
            end

            # do permutations
            if found_gene
                if var.class != 'D' && var.class != 'I'
                    substitution = DNASequence(DNASequence(string(var.variant)))
                    mut_gene_pos = mut_pos - gene_sidx[gene_idx] + 1
                    perm_seq[mut_gene_pos] = substitution
                else
                    if verbose
                        println(STDOUT, "INFO: TRANSLATE MUTATIONS, indel found: ", var)
                        println(STDOUT, "INFO: TRANSLATE MUTATIONS, skipping this cluster: $id")
                    end
                    frame_shifts += 1
                    found_indel = true
                    break
                end
            end
        end

        if found_gene && !(found_indel)
            if strands[gene_idx] == "-"
                BioSequences.reverse_complement!(perm_seq)
            end

            mut_rna = convert(RNASequence, perm_seq)

            mut_protein = translate(mut_rna, code=BioSequences.ncbi_trans_table[genetic_code], allow_ambiguous_codons=true)

            # Comparing WT vs mutated protein and extracting
            protein = proteins[gene_idx][1]
            if verbose
                println(STDOUT, "INFO: TRANSLATE MUTATIONS, mutated sequence: $perm_seq")
                println(STDOUT, "INFO: TRANSLATE MUTATIONS, wt sequence:      ", convert(String, genes_seqs[gene_idx[1]]))
                println(STDOUT, "INFO: TRANSLATE MUTATIONS, mutated protein sequence: $mut_protein")
                println(STDOUT, "INFO: TRANSLATE MUTATIONS, wt protein sequence:      $protein")
            end
            aa_subs = compare_seq(protein, mut_protein; verbose=verbose)

            if haskey(aa_substitutions, chr)
                push!(aa_substitutions[chr], aa_subs)
            else
                aa_substitutions[chr] = [aa_subs]
            end
        end
    end

    return aa_substitutions, frame_shifts
end


function compare_seq(seqa::Any, seqb::Any; verbose::Bool=false)

    muts = String[]
    ix = 0
    for i in seqa
        ix += 1
        if i == seqb[ix]
            continue
        else
            mut = string(i)*string(ix)*string(seqb[ix])
            push!(muts, mut)
        end
    end

    if verbose
        println(STDOUT, "INFO: COMPARE SEQ, found translated mutations: $muts")
    end

    return muts
end


function mapp_mutations(vars::Dict{String, MutationsCluster}; cluster::Bool=true)

    if cluster
        mapp_count = Dict{Array{Mutation,1},Int64}()
    else
        mapp_count = Dict{Mutation,Int64}()
    end

    for id in keys(vars)
        mutations = collect(keys(vars[id].mutation_count))
        if cluster
            if mutations in keys(mapp_count)
                mapp_count[mutations] += 1
            else
                mapp_count[mutations] = 1
            end
        else
            for mutation in mutations
                if mutation in keys(mapp_count)
                    mapp_count[mutation] += 1
                else
                    mapp_count[mutation] = 1
                end
            end
        end
    end

    return mapp_count
end


function parse_chr(id::String)

    return String(split(id,":")[1])
end


function stringify(mutation::Mutation)

    mutation_str = string(mutation.reference)
    mutation_str = mutation_str * string(mutation.mStart)
    mutation_str = mutation_str * string(mutation.variant)

    return mutation_str
end


function conf_coverage!(Coverage::Coverage, vars::Dict{String,MutationsCluster})

    for id in keys(vars)
        chr = parse_chr(id)
        start = vars[id].c_start
        end_ = vars[id].c_start + vars[id].R1_length - 1
        try
            Coverage.coverage[chr][start:end_] += 1
        catch error
            println(STDERR, "ERROR: $error")
            println(STDERR, "ERROR: ", vars[id])
        end
        end_ = vars[id].c_end
        start = vars[id].c_end - vars[id].R2_length + 1
        try
            Coverage.coverage[chr][start:end_] += 1
        catch error
            println(STDERR, "ERROR: $error")
            println(STDERR, "ERROR: ", vars[id])
        end
    end

    return Coverage
end


function transvers(vars::Dict{String,MutationsCluster})

    types = Dict(
        "AG"=>0,
        "GA"=>0,
        "TC"=>0,
        "CT"=>0,
        "AC"=>0,
        "CA"=>0,
        "GT"=>0,
        "TG"=>0,
        "AT"=>0,
        "TA"=>0,
        "CG"=>0,
        "GC"=>0
        )
    trans = 0
    vers = 0
    for id in keys(vars)
        for var in keys(vars[id].mutation_count)
            if var.class == 'M'
                if var.reference == "A" && var.variant == "G"
                    trans += 1
                    types["AG"] += 1
                elseif var.reference == "G" && var.variant == "A"
                    trans += 1
                    types["GA"] += 1
                elseif var.reference == "T" && var.variant == "C"
                    trans += 1
                    types["TC"] += 1
                elseif var.reference == "C" && var.variant == "T"
                    trans += 1
                    types["CT"] += 1
                elseif var.reference == "A" && var.variant == "C"
                    vers += 1
                    types["AC"] += 1
                elseif var.reference == "C" && var.variant == "A"
                    vers += 1
                    types["CA"] += 1
                elseif var.reference == "G" && var.variant == "T"
                    vers += 1
                    types["GT"] += 1
                elseif var.reference == "T" && var.variant == "G"
                    vers += 1
                    types["TG"] += 1
                elseif var.reference == "A" && var.variant == "T"
                    vers += 1
                    types["AT"] += 1
                elseif var.reference == "T" && var.variant == "A"
                    vers += 1
                    types["TA"] += 1
                elseif var.reference == "C" && var.variant == "G"
                    vers += 1
                    types["CG"] += 1
                elseif var.reference == "G" && var.variant == "C"
                    vers += 1
                    types["GC"] += 1
                end
            end
        end
    end

    return trans, vers, types
end


function mapp_single_aa_subs(arrays::Array{Array{String,1}})

    aa_subs = Dict{String, UInt64}()
    for array in arrays
        for aa_sub in array
            if haskey(aa_subs, aa_sub)
                aa_subs[aa_sub] += 1
            else
                aa_subs[aa_sub] = 1
            end
        end
    end

    return aa_subs
end

end  # module
