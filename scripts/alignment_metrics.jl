#!/usr/bin/env julia

using ArgParse.ArgParseSettings
using ArgParse.@add_arg_table
using ArgParse.parse_args
using BioAlignments
using DataStructures
using CSV
using DataFrames
using StatsBase
using GenomicFeatures

include("./CustomJuliaModules/WriteDataFrame.jl")
using WriteDataFrame

function isinregion(pos::Int64,chr::String,
    exclude_reg::Dict{String,Array{String,1}})::Bool

    if haskey(exclude_reg,chr)
        fr = parse(Int64,split(exclude_reg[chr][1],":")[1])
        lr = parse(Int64,split(exclude_reg[chr][end],":")[2])
        if pos >= fr && pos <= lr
            for intv in exclude_reg[chr]
                pos_ex = split(intv,":")

                if pos > parse(Int64,pos_ex[1]) && pos < parse(Int64,pos_ex[2])
                    return true
                end
            end
        end
    end
    return false
end

function bamcounter(bam::String, win::Int64, bed::Any, is_se::Bool)
    data = DataFrame(Sample=String[], TotalReads=Int64[], PrimerAlignedBases=Int64[],
                     OtherBases=Int64[], IsMapped=Int64[], MutAndInDel=Int64[],
                     IsPrimary=Int64[], MapQ0=Int64[], MapQ1=Int64[], MapQ5=Int64[],
                     MapQ10=Int64[], MapQ20=Int64[], MapQ40=Int64[], MapQup40=Int64[],
                     HardClip=Int64[], SoftClip=Int64[], CIGARXID=Int64[],
                     PairedReads=Int64[], MappedProperPairs=Int64[],
                     UnmappedReads=Int64[], MatesUnmapped=Int64[],
                     PCROrOpticalDuplicates=Int64[], SuplementaryAlignment=Int64[],
                     SecondaryAlignment=Int64[])
    starts_data = DataFrame()

    reader = open(BAM.Reader, bam)
    record = BAM.Record()
    prime_aligned_bases = Int64(0)
    other_bases = Int64(0)
    ct = Int64(0)
    ismampped_ct = Int64(0)
    mut_and_indels_ct = Int64(0)
    isprimarry_ct = Int64(0)
    map_q0_ct = Int64(0)
    map_q1_ct = Int64(0)
    map_q5_ct = Int64(0)
    map_q10_ct = Int64(0)
    map_q20_ct = Int64(0)
    map_q40_ct = Int64(0)
    map_up40_ct = Int64(0)
    cigar_hard_clip_ct = Int64(0)
    cigar_soft_clip_ct = Int64(0)
    cigar_seq_mismatch_ct = Int64(0)
    paired_reads_ct = Int64(0)
    proper_pairs_mapped_ct = Int64(0)
    unmapped_reads_ct = Int64(0)
    mates_unmapped_ct = Int64(0)
    read_PCR_or_optical_duplicate_ct = Int64(0)
    suplementary_alignment_ct = Int64(0)
    secondary_alignment_ct = Int64(0)
    mate1 = false
    mate2 = false
    posr = Int64(0)
    posl = Int64(0)
    totalPairedLen = Int64(0)
    totalOverlap = Int64(0)
    overlapPer = Float64(0)
    treads = 0
    wreads = 0
    m2posr = 0
    m1posl = 0
    mate2tmp = ""
    start_bias = Dict{String,Int64}()
    win_nr = Int64(0)
    bed_f = false
    exclude_reg = Dict{String,Array{String,1}}()

    ## Reads excluded regions bed file, if exist.

    if bed != nothing
        bed_f = true
        reader_bed = open(BED.Reader, bed)

        for br in reader_bed
            chr_n = seqname(br)
            intv = string(leftposition(br)) * ":" * string(rightposition(br))

            if haskey(exclude_reg,chr_n)
                push!(exclude_reg[chr_n],intv)
            else
                exclude_reg[chr_n]=[intv]
            end
        end
        close(reader_bed)
    end

    while !eof(reader)
        read!(reader, record)
        flag = BAM.flag(record)
        cigar = BAM.cigar(record)
        ismapped = BAM.ismapped(record)

        ## Overlap calculations
        if ismapped
            cigar = BAM.cigar(record)
            insdel = contains(cigar, "I") || contains(cigar, "D")

            ## Check for first in pair, primary, paired read
            if flag&900 == 0 && flag&3 == 3 && !insdel
                posr = BAM.rightposition(record)
                m1posl = BAM.position(record)
                name = BAM.tempname(record)
                mate1 = true
                mate1tmp = BAM.tempname(record)
                mate1len = posr - m1posl

                if mate2 && m1posl < posl && mate1tmp == mate2tmp
                    lenpair = mate1len + mate2len
                    overlap = posr - posl

                    if lenpair >= overlap
                        totalPairedLen += lenpair

                        if (overlap > 0)
                            totalOverlap += overlap
                        end
                        treads += 1
                    else
                        wreads += 1
                    end
                    mate1 = false
                    mate2 = false
                    posr = 0
                    posl = 0
                    mate1len = 0
                    mate2len = 0
                    m2posr = 0
                    m1posl = 0
                end
            end

            ## Check for second in pair, primary, paired read
            if flag&836 == 0 && flag&3 == 3 && !insdel
                posl = BAM.position(record)
                m2posr = BAM.rightposition(record)
                mate2 = true
                mate2len = m2posr - posl
                mate2tmp = BAM.tempname(record)

                if mate1 && m1posl < posl && mate1tmp == mate2tmp
                    lenpair = mate1len + mate2len
                    overlap = posr - posl

                    if lenpair >= overlap
                        totalPairedLen += lenpair

                        if (overlap > 0)
                            totalOverlap += overlap
                        end
                        treads += 1
                    else
                        wreads += 1
                    end
                    mate1 = false
                    mate2 = false
                    posr = 0
                    posl = 0
                    mate1len = 0
                    mate2len = 0
                    m2posr = 0
                    m1posl = 0
                end
            end
        end

        if totalPairedLen != 0
            overlapPer = round((totalOverlap / totalPairedLen) * 100,2)
        end

        if flag&1 == 1
            paired_reads_ct += 1
        end

        if flag&3 == 3
            proper_pairs_mapped_ct += 1
        end

        if flag&4 == 4
            unmapped_reads_ct += 1
        end

        if flag&9 == 9
            mates_unmapped_ct += 1
        end

        if flag&1024 == 1024
            read_PCR_or_optical_duplicate_ct += 1
        end

        if flag&2048 == 2048
            suplementary_alignment_ct += 1
        end

        if flag&256 == 256
            secondary_alignment_ct += 1
        end

        if flag&2304 == 0
            isprimarry_ct += 1
        end

        if ismapped
            ismampped_ct += 1
        end

        try md_taq = record["MD"]
            splited_MD = split(replace(record["MD"], r"[\\^]*[ACGT]+", x -> "_"*x*"_"), "_")
            if length(splited_MD) > 1
                mut_and_indels_ct += 1
            end
        end

        try mapq = UInt8((record.bin_mq_nl >> 8) & 0xff)
            if mapq == 0
                map_q0_ct += 1
            end

            if mapq in 0:1
                map_q1_ct += 1
            end

            if mapq in 1:5
                map_q5_ct += 1
            end

            if mapq in 5:10
                map_q10_ct += 1
            end

            if mapq in 10:20
                map_q20_ct += 1
            end

            if mapq in 20:40
                map_q40_ct += 1
            end

            if mapq > 40
                map_up40_ct += 1
            end
        end

        if (flag&4 == 0) && (flag&256 == 0) && (flag&2048 == 0)  # filter out only mapped and primary alignments
            prime_aligned_bases += BAM.seqlength(record)
        else
            other_bases += BAM.seqlength(record)
        end

        ## Start_bais calculation

        if (flag&4 == 0) && (flag&256 == 0) && (flag&2048 == 0)
            chr = BAM.refname(record)
            isNotChrUn = !contains(chr, "chrUn")

            if ((flag&64 == 64 && flag&3 == 3) || is_se) && isNotChrUn
                pos = BAM.position(record)
                in_reg = false

                if bed_f
                    in_reg = isinregion(pos,chr,exclude_reg)
                end

                if !in_reg
                    win_nr = floor(Int64,pos/win)
                    win_nr_k = chr * ":" * string(win_nr)

                    if haskey(start_bias, win_nr_k)
                        start_bias[win_nr_k] += 1
                    else
                        start_bias[win_nr_k] = 1
                    end
                end
            end
        end

        if (search(uppercase(cigar), 'X') > 0 ||
            search(uppercase(cigar), 'I') > 0 ||
            search(uppercase(cigar), 'D') > 0)
            cigar_seq_mismatch_ct += 1
        end

        if search(uppercase(cigar), 'S') > 0
            cigar_soft_clip_ct += 1
        end

        if search(uppercase(cigar), 'H') > 0
            cigar_hard_clip_ct += 1
        end


        ct += 1
        if mod(ct,10000) == 0
            print(STDERR, "PARSED $ct BAM records \r")
        end
    end

    push!(data, [splitext(basename(bam))[1], ct, prime_aligned_bases,
                 other_bases, ismampped_ct, mut_and_indels_ct, isprimarry_ct, map_q0_ct,
                 map_q1_ct, map_q5_ct, map_q10_ct, map_q20_ct, map_q40_ct,
                 map_up40_ct, cigar_hard_clip_ct, cigar_soft_clip_ct,
                 cigar_seq_mismatch_ct, paired_reads_ct, proper_pairs_mapped_ct,
                 unmapped_reads_ct, mates_unmapped_ct, read_PCR_or_optical_duplicate_ct,
                 suplementary_alignment_ct, secondary_alignment_ct])

    if isempty(start_bias)
        starts_data[:NrOfStarts] = 0
        starts_data[:Frequency] = 0
    else
        dataForLogs = DataFrame()
        dataForLogs[:Chr] = collect(keys(start_bias))
        dataForLogs[:NrOfStarts] = collect(values(start_bias))
        logN = replace(bam,".bam" => "_starts_bias.csv")
        CSV.write(logN, dataForLogs, delim=',')
        count_uniq = countmap(collect(values(start_bias)))
        starts_data[:NrOfStarts] = collect(keys(count_uniq))
        starts_data[:Frequency] = collect(values(count_uniq))
        sort!(starts_data)
    end
    return data, starts_data, overlapPer
end


function main(args)
    arg_parse_settings = ArgParseSettings(description="Calculates aligment metrics from bam flags")
    @add_arg_table arg_parse_settings begin
        "--bam", "-b"
            arg_type = String
            help = "Bam file."
            required = true
            dest_name = "bam"

        "--win", "-w"
            arg_type = Int64
            help = "Win size for starts bias"
            required = true
            dest_name = "win"

        "--out", "-o"
            arg_type = String
            help = "directory and file name as one string for aligment metrics."
            required = true
            dest_name = "outfile"

        "--bed", "-e"
            arg_type = String
            help = "exclude regions in bed for start bias calculaions"
            required = false
            dest_name = "bed"

        "--outstart", "-s"
            arg_type = String
            help = "directory and file name as one string for starts bias."
            required = true
            dest_name = "outstart"

        "--overlap-out", "-v"
            arg_type = String
            help = "directory and file name as one string."
            required = true
            dest_name = "overlap-out"

        "--single-end", "-g"
            arg_type = Bool
            help = "set true if single-end analysis is done."
            required = true
            dest_name = "is_se"

    end

    parsed_args = parse_args(arg_parse_settings)
    if parsed_args["bed"] == "None"
        parsed_args["bed"] = nothing
    end
    data, starts_data, overlapPer = bamcounter(parsed_args["bam"],parsed_args["win"], parsed_args["bed"], parsed_args["is_se"])
    wdframe(parsed_args["outfile"], data, ',')
    wdframe(parsed_args["outstart"], starts_data, ',')
    # Writing overlap data
    ovdata = DataFrame(Sample=String[], overlapPercent=Float64[])
    push!(ovdata, [splitext(basename(parsed_args["bam"]))[1], overlapPer])
    CSV.write(parsed_args["overlap-out"], ovdata, delim=',')


end

main(ARGS)
