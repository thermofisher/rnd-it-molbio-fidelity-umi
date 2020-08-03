module TestFidelityUMI

using Base.Test
NEW_PATH=join(split(Base.source_path(),"/")[1:end-3],"/")
push!(LOAD_PATH, NEW_PATH)
using BioAlignments
using AutoHashEquals
using DataStructures
using DataFrames
using BioSequences
using BioAlignments
using FidelityUMI

println(LOAD_PATH)

bam_file1 = NEW_PATH*"/FidelityUMI.jl/tests/1_duplicated_subSortbyName.bam"
bam_file2 = NEW_PATH*"/FidelityUMI.jl/tests/2_duplicated_subSortbyName.bam"
reference_file = NEW_PATH*"/FidelityUMI.jl/tests/umi_fidelity.fasta"
gene_csv = NEW_PATH*"/FidelityUMI.jl/tests/genes.csv"

@testset "FidelityUMI" begin
    referenceSruct = ReferenceSequences(reference_file)
    cvCoverage = Coverage(referenceSruct)
    coverage_counter!(cvCoverage, bam_file1; mappingQuality=1)
    clust_coverages = Coverage(referenceSruct)

    @test String["lambda"] == cvCoverage.chromosomes
    @test Dict("lambda"=>379) == cvCoverage.chrom_lengths
    @test Dict("lambda"=>[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 17.0, 17.0, 17.0, 17.0, 17.0, 17.0, 17.0, 17.0, 17.0, 17.0, 17.0, 17.0, 17.0, 17.0, 17.0, 17.0, 17.0, 17.0, 17.0, 17.0, 17.0, 17.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 18.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
              ) == cvCoverage.coverage

    all_variants = Dict(Any[])
    reader = open(BAM.Reader, bam_file2)
    record = BAM.Record()
    while !eof(reader)
        read!(reader, record)
        rname = BAM.tempname(record)
        flag1 = BAM.flag(record)
        variants = call_variants(record; quality=0, verbose=false)

        if flag1&64 == 64
            rname = "R1:"*rname
        elseif flag1&128 == 128
            rname = "R2:"*rname
        end

        if length(variants) > 0
            if rname in keys(all_variants)
                push!(all_variants[rname], variants)
            else
                all_variants[rname] = [variants]
            end
        end
    end
    close(reader)

    @test [
           FidelityUMI.Mutation[FidelityUMI.Mutation('D', "test", 8, 8, "G", "*"),
           FidelityUMI.Mutation('M', "test", 9, 9, "C", "A"),
           FidelityUMI.Mutation('I', "test", 3, 4, "-", "A")]
    ] == all_variants["R1:M00827:86:000000000-BPT2L:1:1:1:1:insert:r1:4A:deletion:8G:mistmatch:r1:C9A_CCCCCCCCCC"]

    @test [
           FidelityUMI.Mutation[FidelityUMI.Mutation('M', "test", 8, 8, "G", "A"),
           FidelityUMI.Mutation('I', "test", 3, 4, "-", "A")]
    ] == all_variants["R1:M00827:86:000000000-BPT2L:1:1:1:1:insert:r1:4A:mistmatch:r1:G8A_CCCCCCCCCC"]

    @test [
           FidelityUMI.Mutation[FidelityUMI.Mutation('M', "test", 8, 8, "G", "C")]
    ] == all_variants["R1:M00827:86:000000000-BPT2L:1:1:1:1:mistmatch:G8C_GGGGGGGGGG"]

    @test [
           FidelityUMI.Mutation[FidelityUMI.Mutation('M', "test", 4, 4, "T", "A")]
    ] == all_variants["R1:M00827:86:000000000-BPT2L:1:1:1:1:mistmatch:r1:T4A_CCCCCCCCCC"]

    @test [
           FidelityUMI.Mutation[FidelityUMI.Mutation('M', "test", 46, 46, "T", "A")]
    ] == all_variants["R2:M00827:86:000000000-BPT2L:1:1:1:1:mistmatch:r2:T5A_CCCCCCCCCC"]

    @test [
           FidelityUMI.Mutation[FidelityUMI.Mutation('M', "test", 4, 4, "T", "A")]
    ] == all_variants["R1:M00827:86:000000000-BPT2L:1:1:1:1:soft:clip:r1:1G:mistmatch:r1:T4A_CCCCCCCCCC"]

    @test [
           FidelityUMI.Mutation[FidelityUMI.Mutation('M', "test", 8, 8, "G", "C")]
    ] == all_variants["R1:M00827:86:000000000-BPT2L:1:1:1:1:mistmatch:G8C_GGGGGGGGGG"]

    @test [
           FidelityUMI.Mutation[FidelityUMI.Mutation('M', "test", 8, 8, "G", "C")]
    ] == all_variants["R1:M00827:86:000000000-BPT2L:1:1:1:2:mistmatch:G8C_GGGGGGGGGG"]

    @test [
           FidelityUMI.Mutation[FidelityUMI.Mutation('M', "test", 8, 8, "G", "C")]
    ] == all_variants["R1:M00827:86:000000000-BPT2L:1:1:1:3:mistmatch:G8C_GGGGGGGGGG"]

    @test [
           FidelityUMI.Mutation[FidelityUMI.Mutation('M', "test", 8, 8, "G", "C")]
    ] == all_variants["R1:M00827:86:000000000-BPT2L:1:1:1:4:mistmatch:G8C_GGGGGGGGGG"]

    @test [
           FidelityUMI.Mutation[FidelityUMI.Mutation('D', "test", 4, 4, "T", "*"),
           FidelityUMI.Mutation('D', "test", 8, 8, "G", "*"),
           FidelityUMI.Mutation('M', "test", 14, 14, "T", "A")]
    ] == all_variants["R1:M00827:86:000000000-BPT2L:1:1:1:1:deletion:r1:4T:8G:mistmatch:r1:T14A_CCCCCCCCCC"]

    @test [
            FidelityUMI.Mutation[FidelityUMI.Mutation('D', "test", 8, 9, "GC", "*"),
            FidelityUMI.Mutation('M', "test", 14, 14, "T", "A"),
            FidelityUMI.Mutation('I', "test", 3, 4, "-", "A")]
    ] == all_variants["R1:M00827:86:000000000-BPT2L:1:1:1:1:insert:r1:4A:deletion:8GC:mistmatch:r1:T14A_CCCCCCCCCC"]

    @test [
            FidelityUMI.Mutation[FidelityUMI.Mutation('M', "test", 14, 14, "T", "A"),
            FidelityUMI.Mutation('I', "test", 7, 8, "-", "AAAA")]
    ] == all_variants["R1:M00827:86:000000000-BPT2L:1:1:1:1:insert:r1:4AAAA:mismatch:T14A_CCCCCCCCCC"]

    @test [
            FidelityUMI.Mutation[FidelityUMI.Mutation('M', "test", 8, 8, "G", "A"),
            FidelityUMI.Mutation('I', "test", 3, 4, "-", "TA")]
    ] == all_variants["R1:M00827:86:000000000-BPT2L:1:1:1:1:insert:r1:4TA:mistmatch:r1:G8A_CCCCCCCCCC"]

    @test [
            FidelityUMI.Mutation[FidelityUMI.Mutation('M', "test", 4, 4, "T", "A"),
            FidelityUMI.Mutation('I', "test", 10, 11, "-", "G")]
    ] == all_variants["R1:M00827:86:000000000-BPT2L:1:1:1:1:mismatch:r1:T4A:insert:r1:10G_CCCCCCCCCC"]

    @test [
            FidelityUMI.Mutation[FidelityUMI.Mutation('M', "test", 4, 4, "T", "A"),
            FidelityUMI.Mutation('M', "test", 7, 7, "T", "A")]
    ] == all_variants["R1:M00827:86:000000000-BPT2L:1:1:1:1:mistmatch:r1:T4A:T7A_CCCCCCCCCC"]

    all_variants = Dict(Any[])
    reader = open(BAM.Reader, bam_file1)
    record = BAM.Record()
    while !eof(reader)
        read!(reader, record)
        rname = BAM.tempname(record)
        variants = call_variants(record; quality=0, verbose=false)
        flag1 = BAM.flag(record)

        if flag1&64 == 64
            rname = "R1:"*rname
        elseif flag1&128 == 128
            rname = "R2:"*rname
        end

        if length(variants) > 0
            if rname in keys(all_variants)
                push!(all_variants[rname], variants)
            else
                all_variants[rname] = [variants]
            end
        end
    end
    close(reader)

    @test [
            FidelityUMI.Mutation[FidelityUMI.Mutation('M', "lambda", 85, 85, "A", "C"),
            FidelityUMI.Mutation('M', "lambda", 100, 100, "T", "G"),
            FidelityUMI.Mutation('M', "lambda", 119, 119, "C", "G")]
    ] == all_variants["R1:M00827:86:000000000-BPT2L:1:2101:17692:1620_ATTACAGGGGACCATGGCAA"]

    @test [
            FidelityUMI.Mutation[FidelityUMI.Mutation('M', "lambda", 100, 100, "T", "G")]
    ] == all_variants["R2:M00827:86:000000000-BPT2L:1:2101:17692:1620_ATTACAGGGGACCATGGCAA"]

    @test [
            FidelityUMI.Mutation[FidelityUMI.Mutation('M', "lambda", 85, 85, "A", "C"),
            FidelityUMI.Mutation('M', "lambda", 100, 100, "T", "G"),
            FidelityUMI.Mutation('M', "lambda", 119, 119, "C", "G")]
    ] == all_variants["R1:M00827:86:000000000-BPT2L:1:2102:17692:1620_ATTACAGGGGACCATGGCAA"]

    @test [
            FidelityUMI.Mutation[FidelityUMI.Mutation('M', "lambda", 100, 100, "T", "G")]
    ] == all_variants["R2:M00827:86:000000000-BPT2L:1:2102:17692:1620_ATTACAGGGGACCATGGCAA"]

    @test [
            FidelityUMI.Mutation[FidelityUMI.Mutation('M', "lambda", 85, 85, "A", "C"),
            FidelityUMI.Mutation('M', "lambda", 100, 100, "T", "G"),
            FidelityUMI.Mutation('M', "lambda", 119, 119, "C", "G")]
    ] == all_variants["R1:M00827:86:000000000-BPT2L:1:2103:17692:1620_ATTACAGGGGACCATGGCAA"]

    @test [
            FidelityUMI.Mutation[FidelityUMI.Mutation('M', "lambda", 100, 100, "T", "G")]
    ] == all_variants["R2:M00827:86:000000000-BPT2L:1:2103:17692:1620_ATTACAGGGGACCATGGCAA"]

    @test [
            FidelityUMI.Mutation[FidelityUMI.Mutation('M', "lambda", 85, 85, "A", "C"),
            FidelityUMI.Mutation('M', "lambda", 100, 100, "T", "G"),
            FidelityUMI.Mutation('M', "lambda", 119, 119, "C", "G")]
    ] == all_variants["R1:M00827:86:000000000-BPT2L:1:2104:17692:1620_ATTACAGGGGACCATGGCAA"]

    @test [
            FidelityUMI.Mutation[FidelityUMI.Mutation('M', "lambda", 100, 100, "T", "G")]
    ] == all_variants["R2:M00827:86:000000000-BPT2L:1:2104:17692:1620_ATTACAGGGGACCATGGCAA"]

    @test [
            FidelityUMI.Mutation[FidelityUMI.Mutation('M', "lambda", 85, 85, "A", "C"),
            FidelityUMI.Mutation('M', "lambda", 100, 100, "T", "G"),
            FidelityUMI.Mutation('D', "lambda", 114, 115, "TA", "*"),
            FidelityUMI.Mutation('M', "lambda", 119, 119, "C", "G"),
            FidelityUMI.Mutation('D', "lambda", 142, 142, "C", "*")]
    ] == all_variants["R1:M00827:86:000000000-BPT2L:1:2105:17692:1620_ATTACAGGGGACCATGGCAA"]

    @test [
            FidelityUMI.Mutation[FidelityUMI.Mutation('D', "lambda", 71, 71, "G", "*"),
            FidelityUMI.Mutation('D', "lambda", 98, 99, "AT", "*"),
            FidelityUMI.Mutation('M', "lambda", 100, 100, "T", "G")]
    ] == all_variants["R2:M00827:86:000000000-BPT2L:1:2105:17692:1620_ATTACAGGGGACCATGGCAA"]

    @test [
            FidelityUMI.Mutation[FidelityUMI.Mutation('D', "lambda", 57, 58, "GA", "*"),
            FidelityUMI.Mutation('M', "lambda", 85, 85, "A", "C"),
            FidelityUMI.Mutation('M', "lambda", 100, 100, "T", "G"),
            FidelityUMI.Mutation('D', "lambda", 114, 115, "TA", "*"),
            FidelityUMI.Mutation('M', "lambda", 119, 119, "C", "G"),
            FidelityUMI.Mutation('D', "lambda", 142, 142, "C", "*")]
    ] == all_variants["R1:M00827:86:000000000-BPT2L:1:2106:17692:1620_ATTACAGGGGACCATGGCAA"]

    @test [
            FidelityUMI.Mutation[FidelityUMI.Mutation('D', "lambda", 71, 71, "G", "*"),
            FidelityUMI.Mutation('D', "lambda", 98, 99, "AT", "*"),
            FidelityUMI.Mutation('M', "lambda", 100, 100, "T", "G"),
            FidelityUMI.Mutation('D', "lambda", 143, 144, "AA", "*")]
    ] == all_variants["R2:M00827:86:000000000-BPT2L:1:2106:17692:1620_ATTACAGGGGACCATGGCAA"]

    @test [
            FidelityUMI.Mutation[FidelityUMI.Mutation('D', "lambda", 57, 58, "GA", "*"),
            FidelityUMI.Mutation('M', "lambda", 85, 85, "A", "C"),
            FidelityUMI.Mutation('M', "lambda", 100, 100, "T", "G"),
            FidelityUMI.Mutation('D', "lambda", 114, 115, "TA", "*"),
            FidelityUMI.Mutation('M', "lambda", 119, 119, "C", "G"),
            FidelityUMI.Mutation('D', "lambda", 142, 142, "C", "*"),
            FidelityUMI.Mutation('I', "lambda", 45, 46, "-", "AA")]
    ] == all_variants["R1:M00827:86:000000000-BPT2L:1:2107:17692:1620_ATTACAGGGGACCATGGCAA"]

    @test [
            FidelityUMI.Mutation[FidelityUMI.Mutation('D', "lambda", 71, 71, "G", "*"),
            FidelityUMI.Mutation('D', "lambda", 98, 99, "AT", "*"),
            FidelityUMI.Mutation('M', "lambda", 100, 100, "T", "G")]
    ] == all_variants["R2:M00827:86:000000000-BPT2L:1:2107:17692:1620_ATTACAGGGGACCATGGCAA"]

    @test [
            FidelityUMI.Mutation[FidelityUMI.Mutation('M', "lambda", 100, 100, "T", "G"),
            FidelityUMI.Mutation('D', "lambda", 114, 115, "TA", "*"),
            FidelityUMI.Mutation('M', "lambda", 119, 119, "C", "G"),
            FidelityUMI.Mutation('D', "lambda", 142, 142, "C", "*")]
    ] == all_variants["R1:M00827:86:000000000-BPT2L:1:2108:17692:1620_ATTACAGGGGACCATGGCAA"]

    all_variants = Dict(Any[])
    reader = open(BAM.Reader, bam_file1)
    record = BAM.Record()
    while !eof(reader)
        read!(reader, record)
        rname = BAM.tempname(record)
        variants = call_variants(record; quality=1, verbose=false)
        flag1 = BAM.flag(record)

        if flag1&64 == 64
            rname = "R1:"*rname
        elseif flag1&128 == 128
            rname = "R2:"*rname
        end

        if length(variants) > 0
            if rname in keys(all_variants)
                push!(all_variants[rname], variants)
            else
                all_variants[rname] = [variants]
            end
        end
    end
    close(reader)

    @test [
            FidelityUMI.Mutation[FidelityUMI.Mutation('M', "lambda", 100, 100, "T", "G")]
    ] == all_variants["R1:M00827:86:000000000-BPT2L:1:2109:17692:1620_ATTACAGGGGACCATGGCAA"]

    @test true == check_operations([BioAlignments.Operation(OP_MATCH)])
    @test true == check_operations([BioAlignments.Operation(OP_INSERT)])
    @test true == check_operations([BioAlignments.Operation(OP_DELETE)])
    @test true == check_operations([BioAlignments.Operation(OP_SOFT_CLIP)])
    @test true == check_operations([BioAlignments.Operation(OP_HARD_CLIP)])
    @test false == check_operations([BioAlignments.Operation(OP_SKIP)])
    @test false == check_operations([BioAlignments.Operation(OP_PAD)])
    @test false == check_operations([BioAlignments.Operation(OP_SEQ_MATCH)])
    @test false == check_operations([BioAlignments.Operation(OP_SEQ_MISMATCH)])
    @test false == check_operations([BioAlignments.Operation(OP_BACK)])
    @test false == check_operations([BioAlignments.Operation(OP_START)])


    (errors,
    mut_coverages,
    total_bases1,
    total_reads,
    filtered_reads,
    filtered_R1,
    filtered_R2,
    clusts_with_vars
    ) = collect_variants(bam_file1, clust_coverages, 0, 10, 0, "_", 2, false, "BX"; verbose=false)

    @test String["lambda"] == clust_coverages.chromosomes
    @test Dict("lambda"=>379) == clust_coverages.chrom_lengths
    @test Dict("lambda"=>[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                ) == clust_coverages.coverage

    (true_errors,
    err_distr,
    clust_size_distr,
    found_mutations,
    found_insertions,
    found_deletions,
    found_ns,
    total_clusters,
    conf_bases,
    clusts_with_vars,
    passing_vars,
    cluster_total_bp,
    conf_coverage
    ) = filter_variants(errors, mut_coverages, cvCoverage, referenceSruct,
                        2, 0.9; verbose=false)

    test_conf_coverage = Dict("lambda"=>UInt32[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    )

    @test test_conf_coverage == conf_coverage.coverage

    test_errors = Dict(
                "lambda:33:164:ATTACAGGGGACCATGGCAA"=>FidelityUMI.MutationsCluster(
                    Dict(
                        FidelityUMI.Mutation('D', "lambda", 71, 71, "G", "*")=>3,
                        FidelityUMI.Mutation('D', "lambda", 143, 144, "AA", "*")=>1,
                        FidelityUMI.Mutation('D', "lambda", 57, 58, "GA", "*")=>2,
                        FidelityUMI.Mutation('M', "lambda", 100, 100, "T", "G")=>3,
                        FidelityUMI.Mutation('D', "lambda", 114, 115, "TA", "*")=>3,
                        FidelityUMI.Mutation('D', "lambda", 142, 142, "C", "*")=>3,
                        FidelityUMI.Mutation('M', "lambda", 85, 85, "A", "C")=>3,
                        FidelityUMI.Mutation('M', "lambda", 119, 119, "C", "G")=>3,
                        FidelityUMI.Mutation('D', "lambda", 98, 99, "AT", "*")=>3,
                        FidelityUMI.Mutation('I', "lambda", 45, 46, "-", "AA")=>1
                    ),
                    3,
                    132,
                    92,
                    10,
                    [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0],
                    [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0],
                    33,
                    164
                ),
                "lambda:63:164:ATTACAGGGGACCATGGCAA"=>FidelityUMI.MutationsCluster(
                    Dict(
                        FidelityUMI.Mutation('D', "lambda", 71, 71, "G", "*")=>1,
                        FidelityUMI.Mutation('M', "lambda", 100, 100, "T", "G")=>1,
                        FidelityUMI.Mutation('D', "lambda", 114, 115, "TA", "*")=>1,
                        FidelityUMI.Mutation('D', "lambda", 142, 142, "C", "*")=>1,
                        FidelityUMI.Mutation('D', "lambda", 98, 99, "AT", "*")=>1,
                        FidelityUMI.Mutation('M', "lambda", 119, 119, "C", "G")=>1
                    ),
                    1,
                    80,
                    80,
                    6,
                    [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
                    [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
                    63,
                    164
                ),
                "lambda:33:154:ATTACAGGGGACCATGGCAA"=>FidelityUMI.MutationsCluster(
                    Dict(
                        FidelityUMI.Mutation('M', "lambda", 100, 100, "T", "G")=>5,
                        FidelityUMI.Mutation('M', "lambda", 85, 85, "A", "C")=>5,
                        FidelityUMI.Mutation('M', "lambda", 119, 119, "C", "G")=>5
                    ),
                    5,
                    122,
                    92,
                    3,
                    [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0],
                    [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0],
                    33,
                    154
                )
            )


    @test test_errors == errors
    @test 18 == total_reads
    @test 1890 == total_bases1
    @test 426 == conf_bases
    @test 1706 == cluster_total_bp
    @test Dict("lambda"=>Dict(
        FidelityUMI.Mutation('D', "lambda", 71, 71, "G", "*")=>4,
        FidelityUMI.Mutation('D', "lambda", 143, 144, "AA", "*")=>1,
        FidelityUMI.Mutation('D', "lambda", 57, 58, "GA", "*")=>2,
        FidelityUMI.Mutation('M', "lambda", 100, 100, "T", "G")=>18,
        FidelityUMI.Mutation('D', "lambda", 114, 115, "TA", "*")=>4,
        FidelityUMI.Mutation('D', "lambda", 142, 142, "C", "*")=>4,
        FidelityUMI.Mutation('M', "lambda", 85, 85, "A", "C")=>8,
        FidelityUMI.Mutation('M', "lambda", 119, 119, "C", "G")=>9,
        FidelityUMI.Mutation('D', "lambda", 98, 99, "AT", "*")=>4,
        FidelityUMI.Mutation('I', "lambda", 45, 46, "-", "AA")=>1)
        ) == mut_coverages

    @test 8 == true_errors
    @test Float32[1.0, 0.33, 0.67, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.33, 1.0, 1.0, 1.0] == err_distr
    @test [3, 1, 5] == clust_size_distr
    @test 4 == found_mutations
    @test 0 == found_insertions
    @test 4 == found_deletions
    @test 0 == found_ns
    @test 3 == total_clusters

    test_passing_vars = Dict(
        "lambda:33:164:ATTACAGGGGACCATGGCAA"=>FidelityUMI.MutationsCluster(
            Dict(
                FidelityUMI.Mutation('D', "lambda", 71, 71, "G", "*")=>1,
                FidelityUMI.Mutation('D', "lambda", 114, 115, "TA", "*")=>1,
                FidelityUMI.Mutation('D', "lambda", 142, 142, "C", "*")=>1,
                FidelityUMI.Mutation('M', "lambda", 85, 85, "A", "C")=>1,
                FidelityUMI.Mutation('M', "lambda", 119, 119, "C", "G")=>1,
                FidelityUMI.Mutation('D', "lambda", 98, 99, "AT", "*")=>1
            ),
            3,
            132,
            92,
            6,
            [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0],
            [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0],
            33,
            164
        ),
        "lambda:33:154:ATTACAGGGGACCATGGCAA"=>FidelityUMI.MutationsCluster(
            Dict(
                FidelityUMI.Mutation('M', "lambda", 85, 85, "A", "C")=>1,
                FidelityUMI.Mutation('M', "lambda", 119, 119, "C", "G")=>1
            ),
            5,
            122,
            92,
            2,
            [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0],
            [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0],
            33,
            154
        )
    )

    @test test_passing_vars == passing_vars
    aa_subs, frame_shifts = translate_mutations(referenceSruct, gene_csv, passing_vars; verbose=false)
    @test Dict("lambda"=>Array{String,1}[String["Q27P", "C38W"]]) == aa_subs
    @test 1 == frame_shifts

    test_stringify = [
        "G71*",
        "AA143*",
        "GA57*",
        "T100G",
        "TA114*",
        "C142*",
        "A85C",
        "C119G",
        "AT98*",
        "-45AA",
        "G71*",
        "T100G",
        "TA114*",
        "C142*",
        "AT98*",
        "C119G",
        "T100G",
        "A85C",
        "C119G",
    ]

    ct = 0
    for id in keys(test_errors)
        for var in keys(test_errors[id].mutation_count)
            ct += 1
            @test test_stringify[ct] == stringify(var)
        end
    end

    test_mapp_mutations = Dict(
            FidelityUMI.Mutation[FidelityUMI.Mutation('M', "lambda", 100, 100, "T", "G"),
            FidelityUMI.Mutation('M', "lambda", 85, 85, "A", "C"),
            FidelityUMI.Mutation('M', "lambda", 119, 119, "C", "G")]=>1,
            FidelityUMI.Mutation[FidelityUMI.Mutation('D', "lambda", 71, 71, "G", "*"),
            FidelityUMI.Mutation('D', "lambda", 143, 144, "AA", "*"),
            FidelityUMI.Mutation('D', "lambda", 57, 58, "GA", "*"),
            FidelityUMI.Mutation('M', "lambda", 100, 100, "T", "G"),
            FidelityUMI.Mutation('D', "lambda", 114, 115, "TA", "*"),
            FidelityUMI.Mutation('D', "lambda", 142, 142, "C", "*"),
            FidelityUMI.Mutation('M', "lambda", 85, 85, "A", "C"),
            FidelityUMI.Mutation('M', "lambda", 119, 119, "C", "G"),
            FidelityUMI.Mutation('D', "lambda", 98, 99, "AT", "*"),
            FidelityUMI.Mutation('I', "lambda", 45, 46, "-", "AA")]=>1,
            FidelityUMI.Mutation[FidelityUMI.Mutation('D', "lambda", 71, 71, "G", "*"),
            FidelityUMI.Mutation('M', "lambda", 100, 100, "T", "G"),
            FidelityUMI.Mutation('D', "lambda", 114, 115, "TA", "*"),
            FidelityUMI.Mutation('D', "lambda", 142, 142, "C", "*"),
            FidelityUMI.Mutation('D', "lambda", 98, 99, "AT", "*"),
            FidelityUMI.Mutation('M', "lambda", 119, 119, "C", "G")]=>1
    )
    #TODO: cluster=false case
    @test test_mapp_mutations == mapp_mutations(test_errors)

    @test String["P27Q", "W38C"] == compare_seq(
                                        "MQFAHKVSGKYRGVAKLEGNTKAKVLPVLATFAYADYWRSA",
                                        "MQFAHKVSGKYRGVAKLEGNTKAKVLQVLATFAYADYCRSA"
                                    )

    cvCoverage = Coverage(referenceSruct)
    test_conf_coverage = Dict("lambda"=>
        UInt32[
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ]
    )

    conf_coverage!(cvCoverage, test_errors)
    @test test_conf_coverage == cvCoverage.coverage

    clust_coverages = Coverage(referenceSruct)
    (errors_paired,
    mut_coverages_paired,
    total_bases1_paired,
    total_reads_paired,
    filtered_reads_paired,
    filtered_R1_paired,
    filtered_R2_paired,
    clusts_with_vars_paired
    ) = collect_variants(bam_file1, clust_coverages, 0, 10, 0, "_", 2, false, "BX"; verbose=false, algo="paired")

    test_errors_paired = Dict(
        "lambda"=>Dict(
            FidelityUMI.Mutation('D', "lambda", 71, 71, "G", "*")=>4,
            FidelityUMI.Mutation('D', "lambda", 143, 144, "AA", "*")=>1,
            FidelityUMI.Mutation('D', "lambda", 57, 58, "GA", "*")=>2,
            FidelityUMI.Mutation('M', "lambda", 100, 100, "T", "G")=>18,
            FidelityUMI.Mutation('D', "lambda", 114, 115, "TA", "*")=>4,
            FidelityUMI.Mutation('D', "lambda", 142, 142, "C", "*")=>4,
            FidelityUMI.Mutation('M', "lambda", 85, 85, "A", "C")=>8,
            FidelityUMI.Mutation('M', "lambda", 119, 119, "C", "G")=>9,
            FidelityUMI.Mutation('D', "lambda", 98, 99, "AT", "*")=>4,
            FidelityUMI.Mutation('I', "lambda", 45, 46, "-", "AA")=>1
            )
        )

end
end
