#!/usr/bin/env bash

# ---------------------------------------------------------------------------- #
# Runs snakemake workflow tests, reports only if error occures. -------------- #
# Runs diff for output files. ------------------------------------------------ #
# Running directory should be root snakemake directory. ---------------------- #
# If --all selected, also builds environment. -------------------------------- #
# During build "--name" environment is cleaned, ".snakemake" is removed ------ #
# and conda caches are cleaned. ---------------------------------------------- #
# Reference files are provided in tests/outputs accordingly to config used for #
# analysis. Checking is done by looping over all reference files, finding -----#
# same files in the new output and comparing with linux diff. ---------------- #
# Logs are saved to tmp_tests directory. ------------------------------------- #
# Arguments: ----------------------------------------------------------------- #
# --all - use to run all tests. Runs build. ---------------------------------- #
# --name - New environment name. Needed if --all. Example 'test' ------------- #
# --path - path to conda environment directory. Example:  -------------------- #
# /media/root/ocean/shared/bin/miniconda/envs/test --------------------------- #
# --umi-indexing - use to run indexing before all UMI tests. ----------------- #
# --umi-fidelity - use to run FidelityUMI.jl workflow tests. ----------------- #
# --clean - remove output files when done. ----------------------------------- #
# --threads - snakemake threads. --------------------------------------------- #
# --parallel-snakes - allowed snakemakes to run in PARALLEL. ----------------- #
# --cluster - run every snake pipeline on the sge cluster. ------------------- #
# --unlock - run to unlock snake crashed runs. ------------------------------- #
# --create-envs - executes snakemake create envs only for all configs. ------- #
# --jenkins - mode for jenkins. ---------------------------------------------- #
# --dryrun - runs all snakemake in dryrun mode. ------------------------------ #
# --sge-env - name of sge pe if running on cluster --------------------------- #
# --sge-queue - name of sge queue if running on cluster ---------------------- #
# --letency-wait - if on cluster snakemake --latency-wait in seconds --------- #
# --verbose - be more verbose about progress --------------------------------- #
# --help - print usage. ------------------------------------------------------ #
# ---------------------------------------------------------------------------- #
# diff params used while comparing files: ------------------------------------ #
# w - ignore all spaces ------------------------------------------------------ #
# N - treat absent files as empty -------------------------------------------- #
# Z - ignore trailing spaces ------------------------------------------------- #
# i - ignore case differencies ----------------------------------------------- #
# E - ignore changes due to tab expansion ------------------------------------ #
# ---------------------------------------------------------------------------- #

RUN_ALL=false
RUN_UMI_FIDELITY=false
RUN_CLEAN=false
RUN_UMI_INDEX=false
RUN_ON_CLUSTER=false
UNLOCK=false
PARALLEL=1
CREATE_ENVS=false
COVERAGE=false
DRYRUN=false
VERBOSE=false

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'
YELLOW='\033[1;33m'
THREADS=8


# ---------------------------------------------------------------------------- #
# Parse arguments ------------------------------------------------------------ #
# ---------------------------------------------------------------------------- #
POSITIONAL=()

while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        -h|--help)
        head -n 41 tests/run_tests.sh | tail -n 40
        exit
        ;;
        -a|--all)
        RUN_ALL=true
        shift
        ;;
        -c|--clean)
        RUN_CLEAN=true
        shift
        ;;
        -t|--threads)
        THREADS=$2
        shift
        shift
        ;;
        -uf|--umi-fidelity)
        RUN_UMI_FIDELITY=true
        shift
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
        -cl|--cluster)
        RUN_ON_CLUSTER=true
        shift
        ;;
        --sq|--sge-queue)
        SGE_QUEUE="-q {$2}"
        shift
        shift
        ;;
        --lw|--latency-wait)
        SNAKE_WAIT="--latency-wait {$2}"
        shift
        shift
        ;;
        --se|--sge-env)
        SGE_PE="{$2}"
        shift
        shift
        ;;
        -ps|--parallel-snakes)
        PARALLEL=$2
        shift
        shift
        ;;
        -un|--unlock)
        UNLOCK=true
        shift
        ;;
        -ce|--create-envs)
        CREATE_ENVS=true
        shift
        ;;
        -ui|--umi-indexing)
        RUN_UMI_INDEX=true
        shift
        ;;
        -v|--verbose)
        VERBOSE=true
        shift
        ;;
        -dr|--dryrun)
        DRYRUN=true
        shift
        ;;
        -cv|--coverage)
        COVERAGE=true
        shift
        ;;
        *)    # unknown option
        POSITIONAL+=("$1") # save it in an array for later
        shift # past argument
        ;;    esac
done

RUN_POLYA=false

echo -e "${GREEN}---- RUNNING PARAMETERS -----"
echo "RUN_ALL          = $RUN_ALL"
echo "RUN_UMI_FIDELITY = $RUN_UMI_FIDELITY"
echo "RUN_CLEAN        = $RUN_CLEAN"
echo "RUN_UMI_INDEX    = $RUN_UMI_INDEX"
echo "THREADS          = $THREADS"
echo "PARALLEL         = $PARALLEL"
echo "UNLOCK           = $UNLOCK"
echo "DRYRUN           = $DRYRUN"
echo "COVERAGE         = $COVERAGE"
echo "RUN_ON_CLUSTER   = $RUN_ON_CLUSTER"
echo -e "--------------------------------${NC}"

if ! [ -z "$POSITIONAL" ]; then
    echo -e "${YELLOW}---- WARNING --------------"
    echo "UNKNOWN PARAMETERS: $POSITIONAL"
    echo -e "--------------------------------${NC}"
fi

# ---------------------------------------------------------------------------- #
# Make tmp directory --------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

if [[ -d tmp_tests ]]; then
    rm -r tmp_tests; mkdir tmp_tests;
else
    mkdir tmp_tests
fi

RULES_DIR="tmp_tests/RULES_COV"
if $COVERAGE; then
    mkdir $RULES_DIR
    grep -Poh 'rule .+:' snakefiles/* | sed 's/rule //' | sed 's/ //g' | sed 's/://' | sort | uniq > $RULES_DIR/ALL_RULES
fi
COUNT_PAR="tmp_tests/COUNT_PARALLEL"
echo 0 > tmp_tests/STARTED_JOBS
echo 0 > tmp_tests/FINISHED_JOBS
echo 0 > tmp_tests/COUNT_PASS
echo 0 > tmp_tests/COUNT_FAIL
echo 0 > tmp_tests/COUNT_TOTAL
echo 0 > tmp_tests/COUNT_MAIN_FAIL
echo 1 > $COUNT_PAR


wait_to_finish(){
    STARTED_JOBS=$(paste -sd+ tmp_tests/STARTED_JOBS | bc)
    FINISHED_JOBS=$(paste -sd+ tmp_tests/FINISHED_JOBS | bc)

    while [ $FINISHED_JOBS -lt $STARTED_JOBS ]; do
        sleep 5
        STARTED_JOBS=$(paste -sd+ tmp_tests/STARTED_JOBS | bc)
        FINISHED_JOBS=$(paste -sd+ tmp_tests/FINISHED_JOBS | bc)
        echo -e -n "\r ${GREEN} STARTED_JOBS: $STARTED_JOBS FINISHED_JOBS: $FINISHED_JOBS ${NC}"
    done
}

wait_parallel(){
    while [ "$(paste -sd+ $COUNT_PAR | bc)" -gt "$PARALLEL" ]; do
        sleep 10
    done
}

if [ -z "$SGE_PE" ]; then
    SGE_PE="smp"
fi
# ---------------------------------------------------------------------------- #
# Test build ----------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

if $RUN_ALL; then

    if [ -z "$name" ]; then
        echo '---------- ERROR ---------------------------'
        echo "--name not provided for environemnt ...."
        echo '--------------------------------------------'
        exit
    fi
    echo -e "${YELLOW} Running conda env build test .... [IN PROGRESS] ${NC}"


    echo 1 >> tmp_tests/COUNT_TOTAL
    if [[ -d .snakemake ]]; then
        echo -e "${RED} INFO: For complete build test .snakemake directory should be removed. ${NC}"
        echo -e "${RED} INFO: Yust make sure that snakemake is not using it at the moment! ${NC}"
        read -p "Do you want to remove '.snakemake' directory [y/n]: " -n 1 -r
        echo    # (optional) move to a new line
        if [[ $REPLY =~ ^[Yy]$ ]]
        then
            echo -e "${YELLOW} WARNING: .snakemake directory will be removed .... [IN PROGRESS] ${NC}"
            rm -r .snakemake;
        fi
    fi

    conda clean --all --yes

    if ! [ -z "$path" ]; then
        if ! bash envs/install_all.sh -n $name -p $path; then
            echo 'ERROR: Failed installation!'
            exit
        fi
    else
        if ! bash envs/install_all.sh -n $name; then
            echo 'ERROR: Failed installation!'
            exit
        fi
    fi
    source activate $name
fi

ENVS_DIR="tmp_tests/CREATE_ENVS"
if  $CREATE_ENVS; then
    mkdir $ENVS_DIR
fi

# wrapper to run the Snakemake.
run_snake(){
    SNAKE_THREADS=$1
    CONFIG_FILE=$2
    PARAMS=$3
    LOG_FILE=$4
    CLUSTER=$5
    MESSAGE=$6
    CONF_NAME=$(basename $CONFIG_FILE .yaml)

    echo 1 >> tmp_tests/COUNT_TOTAL
    if $COVERAGE; then
        echo "snakemake --configfile $CONFIG_FILE --quiet --use-conda -F -n $PARAMS | head -n -1 | tail -n +4 | awk '{print $2}' | sort | uniq >> $RULES_DIR/$CONF_NAME.log"
        snakemake --configfile $CONFIG_FILE  --quiet --use-conda -F -n $PARAMS | head -n -1 | tail -n +4 | awk '{print $2}' | sort | uniq >> $RULES_DIR/$CONF_NAME.log
    elif $CREATE_ENVS; then
        echo "snakemake --configfile $CONFIG_FILE --quiet --use-conda --create-envs-only $PARAMS &> $ENVS_DIR/$CONF_NAME.log"
        snakemake --configfile $CONFIG_FILE --quiet --use-conda --create-envs-only $PARAMS &> $ENVS_DIR/$CONF_NAME.log
    else
        if $CLUSTER; then
            echo -e "${YELLOW} EXECUTING Snakemake: "
            echo -e "${YELLOW} snakemake --quiet --use-conda -j $SNAKE_THREADS -p --configfile $CONFIG_FILE"
            echo -e "${YELLOW} --cluster \"qsub -V -pe $SGE_PE {threads} $SGE_QUEUE -N TEST_{cluster.name} -e {cluster.error} -o {cluster.output} -cwd\""
            echo -e "${YELLOW} --cluster-config cluster.json $PARAMS &> $LOG_FILE ${NC}"

            if snakemake --quiet --use-conda $SNAKE_WAIT -j $SNAKE_THREADS -p -k --configfile $CONFIG_FILE \
                --cluster "qsub -V -pe $SGE_PE {threads} $SGE_QUEUE -N TEST_{cluster.name} -e {cluster.error} -o {cluster.output} -cwd" \
                --cluster-config cluster.json $PARAMS &> $LOG_FILE; then
                echo -e "${GREEN} $MESSAGE .... [OK] ${NC}"
                echo 1 >> tmp_tests/COUNT_PASS
            else
                echo -e "${RED} $MESSAGE .... [ERROR] ${NC}"
                echo 1 >> tmp_tests/COUNT_MAIN_FAIL
                echo 1 >> tmp_tests/COUNT_FAIL
            fi
        else
            echo -e "${YELLOW} EXECUTING Snakemake: "
            echo -e "${YELLOW} snakemake -p -k --configfile $CONFIG_FILE --quiet --use-conda -j $SNAKE_THREADS $PARAMS &> $LOG_FILE ${NC}"

            if snakemake -p -k --configfile $CONFIG_FILE --quiet --use-conda -j $SNAKE_THREADS -k $PARAMS &> $LOG_FILE; then
                echo -e "${GREEN} $MESSAGE .... [OK] ${NC}"
                echo 1 >> tmp_tests/COUNT_PASS
            else
                echo -e "${RED} $MESSAGE .... [ERROR] ${NC}"
                echo 1 >> tmp_tests/COUNT_MAIN_FAIL
                echo 1 >> tmp_tests/COUNT_FAIL
            fi
        fi
    fi
    echo 1 >> tmp_tests/FINISHED_JOBS
}


run_parallel(){
    job=$1
    params=$2
    echo 1 >> tmp_tests/STARTED_JOBS
    wait_parallel
    ($job $params && sed -i '$ d' $COUNT_PAR) & echo 1 >> $COUNT_PAR
}


if $UNLOCK; then
    run_unlock(){
        echo -e "${YELLOW} EXECUTING Snakemake: "
        echo -e "${YELLOW} snakemake --quiet --configfile $1 --unlock ${NC}"
        snakemake --quiet --configfile $1 --unlock
        echo 1 >> tmp_tests/FINISHED_JOBS
    }

    for i in tests/configs/*yaml; do
        run_parallel run_unlock "$i"
    done
fi

if $DRYRUN; then
    for i in tests/configs/*yaml; do
        echo 1 >> tmp_tests/STARTED_JOBS
        fname=$(basename $i)
        run_snake 1 $i "--dryrun" tmp_tests/DRYRUN_${fname:0:-5}.log false "Dry Run: $i"
    done
fi


# ---------------------------------------------------------------------------- #
# Prepare references ----------- --------------------------------------------- #
# ---------------------------------------------------------------------------- #

INDEX_DIR="tmp_tests/INDEXING"
mkdir $INDEX_DIR

# for all UMI

if $RUN_UMI_INDEX || $RUN_ALL; then

    echo 1 >> tmp_tests/STARTED_JOBS
    wait_parallel
    AD_P="get_index"
    (run_snake $THREADS tests/configs/umi_fidelity.yaml "$AD_P" $INDEX_DIR/snakemake_umi_fidelity_index.log \
        $RUN_ON_CLUSTER "umi_fidelity indexing" && sed -i '$ d' $COUNT_PAR) & echo 1 >> $COUNT_PAR

    echo 1 >> tmp_tests/STARTED_JOBS
    wait_parallel
    AD_P="get_index"
    (run_snake $THREADS tests/configs/umi_fidelity2.yaml "$AD_P" $INDEX_DIR/snakemake_umi_fidelity2_index.log \
        $RUN_ON_CLUSTER "umi_fidelity2 indexing" && sed -i '$ d' $COUNT_PAR) & echo 1 >> $COUNT_PAR

    echo 1 >> tmp_tests/STARTED_JOBS
    wait_parallel
    AD_P="get_index"
    (run_snake $THREADS tests/configs/umi_fidelity_paired.yaml "$AD_P" $INDEX_DIR/snakemake_umi_fidelity_paired_index.log \
        $RUN_ON_CLUSTER "umi_fidelity_paired indexing" && sed -i '$ d' $COUNT_PAR) & echo 1 >> $COUNT_PAR
fi

wait_to_finish

# ---------------------------------------------------------------------------- #
# Run umi fidelity small dateset tests --------------------------------------- #
# ---------------------------------------------------------------------------- #

if $RUN_UMI_FIDELITY || $RUN_ALL; then

    UMI_FID="tmp_tests/UMI_FID"
    mkdir $UMI_FID

    run_umi_fidelity1(){
        $VERBOSE && echo -e "${YELLOW} Running tests with umi-fidelity-small dataset ..... [IN PROGRESS] ${NC}"
        $VERBOSE && echo -e "${YELLOW} Running snakemake ........................ [IN PROGRESS] ${NC}"
        AD_P="umi_fidelity"
        run_snake $THREADS tests/configs/umi_fidelity.yaml "$AD_P" $UMI_FID/snakemake_umi_fidelity.log \
            $RUN_ON_CLUSTER "DNA umi dataset umi fidelity workflow"

        $VERBOSE && echo -e "${YELLOW} Checking output with test reference ..... [IN PROGRESS] ${NC}"
        echo 1 >> tmp_tests/STARTED_JOBS
        for i in $(find tests/outputs/umi_fidelity/ \( -iname "*" ! -iname "*.png" ! -iname "*.html" \
                    ! -iname "*.zip" ! -iname "*multiqc*" ! -iname "*config*.yaml" \
                    ! -iname "*.svg"  ! -iname "*.pdf" ! -iname "*.js" ! -iname "*.css" \
                    ! -iname "*.gif" ! -regex '.*/\..*' \) -type f); do
            echo 1 >> tmp_tests/COUNT_TOTAL

            if ( diff -wNZiE "$i" "$(find output_test_umi_fidelity/ -path "*$(echo $i | sed 's/tests\/outputs\/umi_fidelity//')")" ) > $UMI_FID/umi_fidelity_$(basename $i) 2>&1; then
                $VERBOSE && echo -e "${GREEN} $(basename $i) ..... [OK] ${NC}" && echo 1 >> tmp_tests/COUNT_PASS
            else
                $VERBOSE && echo -e "${RED} $(basename $i) ..... [ERROR] ${NC}" && echo 1 >> tmp_tests/COUNT_FAIL
            fi
        done
        echo 1 >> tmp_tests/FINISHED_JOBS
    }

    run_umi_fidelity2(){
        $VERBOSE && echo -e "${YELLOW} Running tests with umi-fidelity-small second dataset ..... [IN PROGRESS] ${NC}"
        $VERBOSE && echo -e "${YELLOW} Running snakemake ........................ [IN PROGRESS] ${NC}"
        AD_P="umi_fidelity"
        run_snake $THREADS tests/configs/umi_fidelity2.yaml "$AD_P" $UMI_FID/snakemake_umi_fidelity2.log \
            $RUN_ON_CLUSTER "DNA umi dataset umi fidelity workflow, second reference"

        $VERBOSE && echo -e "${YELLOW} Checking output with test reference ..... [IN PROGRESS] ${NC}"
        echo 1 >> tmp_tests/STARTED_JOBS
        for i in $(find tests/outputs/umi_fidelity2/ \( -iname "*" ! -iname "*.png" ! -iname "*.html" \
                    ! -iname "*.zip" ! -iname "*multiqc*" ! -iname "*config*.yaml" \
                    ! -iname "*.svg"  ! -iname "*.pdf" ! -iname "*.js" ! -iname "*.css" \
                    ! -iname "*.gif" ! -regex '.*/\..*' \) -type f); do
            echo 1 >> tmp_tests/COUNT_TOTAL

            if ( diff -wNZiE "$i" "$(find output_test_umi_fidelity2/ -path "*$(echo $i | sed 's/tests\/outputs\/umi_fidelity2//')")" ) > $UMI_FID/umi_fidelity2_$(basename $i) 2>&1; then
                $VERBOSE && echo -e "${GREEN} $(basename $i) ..... [OK] ${NC}" && echo 1 >> tmp_tests/COUNT_PASS
            else
                $VERBOSE && echo -e "${RED} $(basename $i) ..... [ERROR] ${NC}" && echo 1 >> tmp_tests/COUNT_FAIL
            fi
        done
        echo 1 >> tmp_tests/FINISHED_JOBS
    }

    run_umi_fidelity3(){
        $VERBOSE && echo -e "${YELLOW} Running tests with umi-fidelity-small paired dataset ..... [IN PROGRESS] ${NC}"
        $VERBOSE && echo -e "${YELLOW} Running snakemake ........................ [IN PROGRESS] ${NC}"

        AD_P="umi_fidelity"
        run_snake $THREADS tests/configs/umi_fidelity_paired.yaml "$AD_P" $UMI_FID/snakemake_umi_fidelity_paired.log \
            $RUN_ON_CLUSTER "DNA umi dataset umi fidelity paired workflow"

        $VERBOSE && echo -e "${YELLOW} Checking output with test reference ..... [IN PROGRESS] ${NC}"
        echo 1 >> tmp_tests/STARTED_JOBS
        for i in $(find tests/outputs/umi_fidelity_paired/ \( -iname "*" ! -iname "*.png" ! -iname "*.html" \
                    ! -iname "*.zip" ! -iname "*multiqc*" ! -iname "*config*.yaml" \
                    ! -iname "*.svg"  ! -iname "*.pdf" ! -iname "*.js" ! -iname "*.css" \
                    ! -iname "*.gif" ! -regex '.*/\..*' \) -type f); do
            echo 1 >> tmp_tests/COUNT_TOTAL

            if ( diff -wNZiE "$i" "$(find output_test_umi_fidelity_paired/ -path "*$(echo $i | sed 's/tests\/outputs\/umi_fidelity_paired//')")" ) > $UMI_FID/umi_fidelity_paired_$(basename $i) 2>&1; then
                $VERBOSE && echo -e "${GREEN} $(basename $i) ..... [OK] ${NC}" && echo 1 >> tmp_tests/COUNT_PASS
            else
                $VERBOSE && echo -e "${RED} $(basename $i) ..... [ERROR] ${NC}" && echo 1 >> tmp_tests/COUNT_FAIL
            fi
        done
        echo 1 >> tmp_tests/FINISHED_JOBS
    }

    run_umi_fidelity4(){
        $VERBOSE && echo -e "${YELLOW} Running FidelityUMI julia unit tests ..... [IN PROGRESS] ${NC}"

        echo 1 >> tmp_tests/COUNT_TOTAL
        if ( julia scripts/FidelityUMI.jl/tests/tests.jl ) > $UMI_FID/FidelityUMI_unit_tests.log 2>&1; then
            $VERBOSE && echo -e "${GREEN} FidelityUMI Unit tests ..... [OK] ${NC}" && echo 1 >> tmp_tests/COUNT_PASS
        else
            $VERBOSE && echo -e "${RED} FidelityUMI Unit tests ..... [ERROR] ${NC}" && echo 1 >> tmp_tests/COUNT_FAIL
        fi
        echo 1 >> tmp_tests/FINISHED_JOBS
    }

    run_umi_fidelity5(){
        $VERBOSE && echo -e "${YELLOW} Running tests with umi-fidelity-small paired bbduk dataset ..... [IN PROGRESS] ${NC}"
        $VERBOSE && echo -e "${YELLOW} Running snakemake ........................ [IN PROGRESS] ${NC}"

        AD_P="umi_fidelity"
        run_snake $THREADS tests/configs/umi_fidelity_paired_bbduk.yaml "$AD_P" $UMI_FID/snakemake_umi_fidelity_paired_bbduk.log \
            $RUN_ON_CLUSTER "DNA umi dataset umi fidelity paired bbduk workflow"

        $VERBOSE && echo -e "${YELLOW} Checking output with test reference ..... [IN PROGRESS] ${NC}"
        echo 1 >> tmp_tests/STARTED_JOBS
        for i in $(find tests/outputs/umi_fidelity_paired/ \( -iname "*" ! -iname "*.png" ! -iname "*.html" \
                    ! -iname "*.zip" ! -iname "*multiqc*" ! -iname "*config*.yaml" \
                    ! -iname "*.svg"  ! -iname "*.pdf" ! -iname "*.js" ! -iname "*.css" \
                    ! -iname "*.gif" ! -regex '.*/\..*' \) -type f); do
            echo 1 >> tmp_tests/COUNT_TOTAL

            if ( diff -wNZiE "$i" "$(find output_test_umi_fidelity_paired_bbduk/ -path "*$(echo $i | sed 's/tests\/outputs\/umi_fidelity_paired//')")" ) > $UMI_FID/umi_fidelity_paired_bbduk_$(basename $i) 2>&1; then
                $VERBOSE && echo -e "${GREEN} $(basename $i) ..... [OK] ${NC}" && echo 1 >> tmp_tests/COUNT_PASS
            else
                $VERBOSE && echo -e "${RED} $(basename $i) ..... [ERROR] ${NC}" && echo 1 >> tmp_tests/COUNT_FAIL
            fi
        done
        echo 1 >> tmp_tests/FINISHED_JOBS
    }

    run_parallel run_umi_fidelity1
    run_parallel run_umi_fidelity2
    run_parallel run_umi_fidelity3
    run_parallel run_umi_fidelity4
    run_parallel run_umi_fidelity5
fi

wait_to_finish

# ---------------------------------------------------------------------------- #
# If env created - remove  --------------------------------------------------- #
# ---------------------------------------------------------------------------- #

if $RUN_ALL; then
    conda env remove -n $name --yes
fi


# ---------------------------------------------------------------------------- #
# Print result  -------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

COUNT_PASS=$(paste -sd+ tmp_tests/COUNT_PASS | bc)
COUNT_FAIL=$(paste -sd+ tmp_tests/COUNT_FAIL | bc)
COUNT_TOTAL=$(paste -sd+ tmp_tests/COUNT_TOTAL | bc)
COUNT_MAIN_FAIL=$(paste -sd+ tmp_tests/COUNT_MAIN_FAIL | bc)

if $COVERAGE; then
    cat $RULES_DIR/*.log | sort | uniq > $RULES_DIR/USED_RULES
    grep -Poh 'rule .+:' snakefiles/* | sed 's/rule //' | sed 's/ //g' | sed 's/://' | sort | uniq > $RULES_DIR/ALL_RULES
    grep -xvf $RULES_DIR/USED_RULES $RULES_DIR/ALL_RULES > $RULES_DIR/LEFT_RULES
else
    if [ $COUNT_FAIL -eq 0 ]; then
        echo
        echo ---------------------------------------------
        echo -e "${GREEN} ALL TESTS PASSED: ${COUNT_PASS} ${NC}";
    else
        echo
        echo ---------------------------------------------
        echo -e "${YELLOW} TOTAL TESTS: ${COUNT_TOTAL}"
        echo -e "${GREEN} TOTAL PASSED: ${COUNT_PASS}"
        echo -e "${RED} TOTAL FAILED: ${COUNT_FAIL} ${NC}"
    fi
    if [ $COUNT_MAIN_FAIL -eq 0 ]; then
        echo -e "${GREEN} ALL MAIN WORKFLOWs PASSED ${NC}"
    else
        echo -e "${RED} MAIN WORKFLOWs FAILED: ${COUNT_MAIN_FAIL} ${NC}"
    fi
fi



# ---------------------------------------------------------------------------- #
# Remove tmp dir ------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
if $RUN_CLEAN && [[ -d tmp_tests ]]; then
    rm -r tmp_tests;
fi

if $RUN_CLEAN && [[ -d output_test_umi_fidelity ]]; then
    rm -r output_test_umi_fidelity;
fi

if $RUN_CLEAN && [[ -d tmp_test_umi_fidelity ]]; then
    rm -r tmp_test_umi_fidelity;
fi

if $RUN_CLEAN && [[ -d output_test_umi_fidelity2 ]]; then
    rm -r output_test_umi_fidelity2;
fi

if $RUN_CLEAN && [[ -d tmp_test_umi_fidelity2 ]]; then
    rm -r tmp_test_umi_fidelity2;
fi

if $RUN_CLEAN && [[ -d output_test_umi_fidelity_paired ]]; then
    rm -r output_test_umi_fidelity_paired;
fi

if $RUN_CLEAN && [[ -d tmp_test_umi_fidelity_paired ]]; then
    rm -r tmp_test_umi_fidelity_paired;
fi
