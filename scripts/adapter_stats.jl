#!/usr/bin/env julia
using GZip
using BioSequences
using ArgParse.ArgParseSettings
using ArgParse.@add_arg_table
using ArgParse.parse_args
using DataStructures
using DataFrames
using CSV

#=
  Calculates adapter-dimer and adapter ratio in total reads.
=#

function adapter_counter(fastq_file::String, csv_out::String, fw_adapt::String)
    adapt_array = []

    f = open(fw_adapt)
    lines = readlines(f)

    for i in lines
        push!(adapt_array, split(i, r"[\s\t]+")[1])
    end

    csv_data::Array{String} = []

    # Adds header.
    data = DataFrame(sample=String[], total_reads=Int64[], adapter=Float64[],
                     adapter_dimer=Int64[])

    # Reads every FASTQ file and calculates its occurence and number of total
    # reads.
    read_counter = 0
    adapter_counter = 0
    adapter_dimer_counter = 0
    stream = GZip.open(fastq_file, "r")

     while !eof(stream)
        name = rstrip(readline(stream),'\n')
        sequence = rstrip(readline(stream),'\n')
        strand = rstrip(readline(stream),'\n')
        quality = rstrip(readline(stream),'\n')
        read_seq = DNASequence(String(sequence))

        approx = 0
        for i in adapt_array
            approx2 = approxsearchindex(read_seq, DNASequence(String(i)), 1)

            if approx < approx2
                approx = approx2
            end
        end

        if approx  > 0
            adapter_counter += 1
        end

        if  approx in 1:4
            adapter_dimer_counter += 1
        end

        read_counter += 1
    end
    close(stream)

    push!(data, [splitext(basename(fastq_file))[1], read_counter,
                 adapter_counter, adapter_dimer_counter])
    CSV.write(csv_out, data)
end


function main(args)
    arg_parse_settings = ArgParseSettings(description="Calculates aligned bases.")
    @add_arg_table arg_parse_settings begin
        "--fastq", "-f"
            arg_type = String
            help = "fastq.gz file."
            required = true
            dest_name = "fastq_f"
        "--out", "-o"
            arg_type = String
            help = "directory and file name as one string."
            required = false
            dest_name = "outfile"
            default = "./output/adapter_metrics_all_samples.csv"
        "--fwd", "-w"
            arg_type = String
            help = "File for adapters."
            required = true
            dest_name = "fwd"
    end

    parsed_args = parse_args(arg_parse_settings)
    adapter_counter(parsed_args["fastq_f"], parsed_args["outfile"],
                    parsed_args["fwd"])

end

main(ARGS)
