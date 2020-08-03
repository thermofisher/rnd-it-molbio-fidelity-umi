#!/usr/bin/env bash
# ---------------------------------------------------------------------------- #
# Run to install julia pkgs on conda environment.  --------------------------- #
# Arguments: ----------------------------------------------------------------- #
# --name | -n - conda env name. Always use since is needed to activate env. -- #
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Parse arguments ------------------------------------------------------------ #
# ---------------------------------------------------------------------------- #
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        -h|--help)
        head -n 10 envs/build_julia.sh | tail 9
        exit
        ;;
        -n|--name)
        name=$2
        shift
        shift
        ;;
        *)    # unknown option
        POSITIONAL+=("$1") # save it in an array for later
        shift # past argument
        ;;    esac
done

echo 'Selected params:' $path $name

if ! [ -z "$POSITIONAL" ]; then
    echo '---------- WARNING --------------'
    echo "UNKNOWN OPTIONS: $POSITIONAL"
    echo '---------------------------------'
fi

if [ -z "$name" ]; then
    echo '---------- ERROR --------------'
    echo "environment --name not provided ...."
    echo '---------------------------------'
    exit
fi

source activate $name
JULIA=julia
$JULIA -e 'Pkg.init()'

if ! $JULIA -e 'using ArgParse' > /dev/null 2>&1; then
    $JULIA -e 'Pkg.add("ArgParse")'
    $JULIA -e 'using ArgParse'
fi

if ! $JULIA -e 'using AutoHashEquals' > /dev/null 2>&1; then
    $JULIA -e 'Pkg.add("AutoHashEquals")'
    $JULIA -e 'using AutoHashEquals'
fi

if ! $JULIA -e 'using BioSequences' > /dev/null 2>&1; then
    $JULIA -e 'Pkg.add("BioSequences", v"0.8.3")'
    $JULIA -e 'using BioSequences'
fi

if ! $JULIA -e 'using BioAlignments' > /dev/null 2>&1; then
    $JULIA -e 'Pkg.add("BioAlignments", v"0.1.0")'
fi

if ! $JULIA -e 'using CSV' > /dev/null 2>&1; then
    $JULIA -e 'Pkg.add("CSV", v"0.1.5")'
    $JULIA -e 'using CSV'
fi

if ! $JULIA -e 'using DataFrames' > /dev/null 2>&1; then
    $JULIA -e 'Pkg.add("DataFrames", v"0.10.1")'
    $JULIA -e 'using DataFrames'
fi

if ! $JULIA -e 'using StatsBase' > /dev/null 2>&1; then
    $JULIA -e 'Pkg.add("StatsBase", v"0.21.0")'
    $JULIA -e 'using StatsBase'
fi

if ! $JULIA -e 'using Lint' > /dev/null 2>&1; then
    $JULIA -e 'Pkg.add("Lint")'
    $JULIA -e 'using Lint'
fi

if ! $JULIA -e 'using LightXML' > /dev/null 2>&1; then
    $JULIA -e 'Pkg.add("LightXML")'
    $JULIA -e 'using LightXML'
fi
