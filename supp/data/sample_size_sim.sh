#!/usr/bin/env bash

set -e

##############################
## Sample Size Simulations ###
##############################


# 60 individuals, 100 generations, 10000 loci, Baseline
SEEDS=(11317 2885 45451 8291 35115 9425 28742 52449 21764 61649 38533)
POP_SIZES="2500 2500 2500"
NUM_POPS="3"
LOCI="10000"
SAMPLES="20 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20"
PREFIX="BASELINE_100GEN_10000LOCI_60D"
SET="baseline"

for i in $(seq 1 10);
do
    if [ ${i} -lt 10 ]; then
        FOLDER=${PREFIX}"_0"${i}
    else
        FOLDER=${PREFIX}"_"${i}
    fi
    mkdir ${FOLDER}
    cd ${FOLDER}
    echo "python ../../scripts/simulate.py -k ${NUM_POPS} -l ${LOCI} -s ${SEEDS[i]} --set ${SET} --sizes ${POP_SIZES} --samples ${SAMPLES}" > README
    python ../../scripts/simulate.py -k ${NUM_POPS} -l ${LOCI} -s ${SEEDS[i]} --set ${SET} --sizes ${POP_SIZES} --samples ${SAMPLES} >> README
    cd ../
done

# 90 individuals, 100 generations, 10000 loci, Baseline
SEEDS=(34695 38071 54193 57405 63422 54504 29971 14959 5767 55561 46705)
POP_SIZES="2500 2500 2500"
NUM_POPS="3"
LOCI="10000"
SAMPLES="30 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 30 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 30"
PREFIX="BASELINE_100GEN_10000LOCI_90D"
SET="baseline"

for i in $(seq 1 10);
do
    if [ ${i} -lt 10 ]; then
        FOLDER=${PREFIX}"_0"${i}
    else
        FOLDER=${PREFIX}"_"${i}
    fi
    mkdir ${FOLDER}
    cd ${FOLDER}
    echo "python ../../scripts/simulate.py -k ${NUM_POPS} -l ${LOCI} -s ${SEEDS[i]} --set ${SET} --sizes ${POP_SIZES} --samples ${SAMPLES}" > README
    python ../../scripts/simulate.py -k ${NUM_POPS} -l ${LOCI} -s ${SEEDS[i]} --set ${SET} --sizes ${POP_SIZES} --samples ${SAMPLES} >> README
    cd ../
done


# 240 individuals, 100 generations, 10000 loci, Baseline
SEEDS=(28248 40671 36886 22927 31432 60353 11198 12169 65093 10316 4431)
POP_SIZES="2500 2500 2500"
NUM_POPS="3"
LOCI="10000"
SAMPLES="80 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 80 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 80"
PREFIX="BASELINE_100GEN_10000LOCI_240D"
SET="baseline"

for i in $(seq 1 10);
do
    if [ ${i} -lt 10 ]; then
        FOLDER=${PREFIX}"_0"${i}
    else
        FOLDER=${PREFIX}"_"${i}
    fi
    mkdir ${FOLDER}
    cd ${FOLDER}
    echo "python ../../scripts/simulate.py -k ${NUM_POPS} -l ${LOCI} -s ${SEEDS[i]} --set ${SET} --sizes ${POP_SIZES} --samples ${SAMPLES}" > README
    python ../../scripts/simulate.py -k ${NUM_POPS} -l ${LOCI} -s ${SEEDS[i]} --set ${SET} --sizes ${POP_SIZES} --samples ${SAMPLES} >> README
    cd ../
done