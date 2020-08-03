#!/usr/bin/env julia

#=
  Reads csv files and concatenates.
  Ex.:
      ./join_csv.jl *.csv output.csv
=#

function main(args)

    csv_data::Array{String} = []
    output_csv::String = args[end]

    # Adds header.
    open(args[1], "r") do f
        lines::Array{String} = readlines(f)
    line_counter::Int = 0
        for line::String in lines
            line_counter += 1
            if line_counter == 1
                        push!(csv_data, line * "\n")
                break
            end

        end
    end

    # reads files and adds lines to csv_data
    for file::String in args[1:end-1]
        # Reads file one by one.
        open(file, "r") do f
            # File parsing and converting to csv.
            lines::Array{String} = readlines(f)
            line_counter::Int = 0

            for line::String in lines
                line_counter += 1

                if line_counter == 1
                    if line * "\n" != csv_data[1]
                        println("Header of input csv files does not match!")
                        exit(1)
                    end
                end

                if line_counter != 1
                    push!(csv_data, line * "\n")
                end
            end
        end
    end

    # Write csv data to a file.
    csv_data_all::String = join(csv_data)
    csv_data_all = replace(csv_data_all, "\t", ",")
    open(output_csv, "w") do f
        write(f, csv_data_all)
    end
end

main(ARGS)
