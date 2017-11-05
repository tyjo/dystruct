#!/usr/bin/env bash

set -e

# Helper script to simulate data. Calls simulate.py with specified parameters

##############################
#### Baseline Simulations ####
##############################

# 120 individuals, 100 generations, 10000 loci, Baseline
SEEDS=(13254 2671 59236 64015 40396 11662 7171 43296 19980 31964 50625)
POP_SIZES="2500 2500 2500"
NUM_POPS="3"
LOCI="10000"
SAMPLES="40 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 40 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 40"
PREFIX="BASELINE_100GEN_10000LOCI_120D"
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

# 120 individuals, 200 generations, 10000 loci, Baseline
SEEDS=(1227 38379 42754 25229 25103 35075 8061 49654 22681 55305 14315)
POP_SIZES="2500 2500 2500"
NUM_POPS="3"
LOCI="10000"
SAMPLES="40 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 40 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 40"
PREFIX="BASELINE_200GEN_10000LOCI_120D"
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

# 120 individuals, 400 generations, 10000 loci, Baseline
SEEDS=(11766 62755 38526 55609 7568 58570 8749 3305 8048 43955 47140)
POP_SIZES="2500 2500 2500"
NUM_POPS="3"
LOCI="10000"
SAMPLES="40 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 40 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 40"
PREFIX="BASELINE_400GEN_10000LOCI_120D"
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

# 120 individuals, 800 generations, 10000 loci, Baseline
SEEDS=(8949 40980 59993 16400 32171 51091 7331 13279 28772 43481 60704)
POP_SIZES="2500 2500 2500"
NUM_POPS="3"
LOCI="10000"
SAMPLES="40 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 40 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 40"
PREFIX="BASELINE_800GEN_10000LOCI_120D"
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





##############################
### Clustered Simulations ####
##############################

# 120 individuals, 100 generations, 10000 loci, Clustered
SEEDS=(28359 18964 6853 46348 27388 24675 25885 14155 50896 42145 24363)
POP_SIZES="2500 2500 2500"
NUM_POPS="3"
LOCI="10000"
SAMPLES="40 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 40 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 40"
PREFIX="CLUSTERED_100GEN_10000LOCI_120D"
SET="clustered"

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

# 120 individuals, 200 generations, 10000 loci, Clustered
SEEDS=(44266 33842 50028 37980 5370 29816 27742 64932 60526 37177 19628)
POP_SIZES="2500 2500 2500"
NUM_POPS="3"
LOCI="10000"
SAMPLES="40 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 40 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 40"
PREFIX="CLUSTERED_200GEN_10000LOCI_120D"
SET="clustered"

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

# 120 individuals, 400 generations, 10000 loci, Clustered
SEEDS=(22641 56926 12472 49215 57401 14251 19083 38355 4127 10631 46723)
POP_SIZES="2500 2500 2500"
NUM_POPS="3"
LOCI="10000"
SAMPLES="40 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 40 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 40"
PREFIX="CLUSTERED_400GEN_10000LOCI_120D"
SET="clustered"

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

# 60 individuals, 800 generations, 10000 loci, Clustered
SEEDS=(60405 48254 28695 29242 40435 40746 19552 63244 39044 1790 26977)
POP_SIZES="2500 2500 2500"
NUM_POPS="3"
LOCI="10000"
SAMPLES="40 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 40 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 40"
PREFIX="CLUSTERED_800GEN_10000LOCI_120D"
SET="clustered"

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