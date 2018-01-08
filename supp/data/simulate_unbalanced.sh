#!/usr/bin/env bash

set -e

# Helper script to simulate data. Calls simulate.py with specified parameters

##############################
#### Unbalanced Simulations ####
##############################

# 23 ancient, 485 modern, 100 generations, 10000 loci, Unbalanced
SEEDS=(20183 5654 24026 2119 818 22383 5538 21753 7275 16921 28609)
POP_SIZES="2500 2500 2500 2500"
NUM_POPS="4"
LOCI="10000"
SAMPLES="1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 11 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 11 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 485"
PREFIX="UNBALANCED_100GEN_10000LOCI_508D"
SET="unbalanced"

for i in $(seq 1 11);
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


# 23 ancient, 485 modern, 200 generations, 10000 loci, Unbalanced
SEEDS=(16294 1327 6038 18967 16895 1337 14125 3348 13859 25469 4778)
POP_SIZES="2500 2500 2500 2500"
NUM_POPS="4"
LOCI="10000"
SAMPLES="1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 11 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 11 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 485"
PREFIX="UNBALANCED_200GEN_10000LOCI_508D"
SET="unbalanced"

for i in $(seq 1 11);
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


# 23 ancient, 485 modern, 400 generations, 10000 loci, Unbalanced
SEEDS=(10173 20409 32536 13016 16512 2768 2807 1983 9913 11382 22601)
POP_SIZES="2500 2500 2500 2500"
NUM_POPS="4"
LOCI="10000"
SAMPLES="1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 11 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 11 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 485"
PREFIX="UNBALANCED_400GEN_10000LOCI_508D"
SET="unbalanced"

for i in $(seq 1 11);
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



# 23 ancient, 485 modern, 800 generations, 10000 loci, Unbalanced
SEEDS=(158 20736 12570 23461 10151 29753 17570 22779 15741 32601 8267)
POP_SIZES="2500 2500 2500 2500"
NUM_POPS="4"
LOCI="10000"
SAMPLES="1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 11 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 11 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 485"
PREFIX="UNBALANCED_800GEN_10000LOCI_508D"
SET="unbalanced"

for i in $(seq 1 11);
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