#!/usr/bin/env bash

# ---------------------------------------------------------------------------- #
# Run to install conda environment. Conda should be preinstalled. ------------ #
# Arguments: ----------------------------------------------------------------- #
# --name | -n - conda env name. Always use since is needed to activate env. -- #
# --path | -p - conda env path. ---------------------------------------------- #
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Parse arguments ------------------------------------------------------------ #
# ---------------------------------------------------------------------------- #
POSITIONAL=()

while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        -h|--help)
        head -n 12 envs/install_all.sh | tail -n 11
        exit
        ;;
        -n|--name)
        name=$2
        shift
        shift
        ;;
        -p|--path)
        path=$2
        shift
        shift
        ;;
        *)    # unknown option
        POSITIONAL+=("$1") # save it in an array for later
        shift # past argument
        ;;    esac
done

echo 'INFO: Selected params: PATH:' $path 'NAME:' $name


if ! [ -z "$POSITIONAL" ]; then
    echo "WARNING: UNKNOWN OPTIONS: $POSITIONAL"
fi

if [ -z "$name" ]; then
    echo "ERROR: --name not provided ...."
    exit
fi


echo 'INFO: INSTALL LOG will be saved to install.log'

if ! [ -z "$path" ]; then
    echo "INFO: Installing conda environment to: $path"
    if ! conda env create -p $path -f envs/environment.yaml > install.log 2>&1; then
        echo 'ERROR: Failed ""conda env create -p "' $path '-f envs/environment.yaml"'
        echo 'INFO: Trying to rewrite ....'
        if ! conda env remove -n $name 2>&1; then
           echo 'ERROR: Failed to remove' $name
           exit
        else
            if ! conda env create -p $path -f envs/environment.yaml >> install.log 2>&1; then
                echo 'ERROR: Failed second time ""conda env create -p "' $path '-f envs/environment.yaml"'
                exit
            fi
        fi
    else
        echo 'SUCCESS: create ENV ' $path
    fi
else
    echo 'INFO: Creating conda environment:' $name '....'
    if ! conda env create -n $name -f envs/environment.yaml > install.log 2>&1; then
        echo 'ERROR: Failed "conda env create -n ' $name '-f envs/environment.yaml"'
        echo 'INFO: Trying to rewrite ....'
        if ! conda env remove -n $name 2>&1; then
           echo 'ERROR: Failed to remove' $name
           exit
        else
            if ! conda env create -n $name -f envs/environment.yaml >> install.log2>&1; then
                echo 'ERROR: Failed second time ""conda env create -n "' $name '-f envs/environment.yaml"'
                exit
            fi
        fi
    else
        echo 'SUCCESS: create ENV ' $name
    fi
fi

if ! source activate $name >> install.log 2>&1; then
    echo 'ERROR: Failed to activate or create environment'
    exit
fi

echo 'INFO: Installing Julia additional packages ....'
if ! bash envs/build_julia_pkgs.sh -n $name >> install.log 2>&1; then
    echo 'ERROR: Failed build julia pakages ....'
    exit
else
    echo 'SUCCESS: Finished build julia packages ....'
fi

echo 'SUCCESS: Finished environemt install....'
