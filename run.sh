#!/usr/bin/env bash
export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/Users/tyjo/dystruct/boost_1_62_0/stage

./bin/dystruct --input ./supp/data/BASELINE_100GEN_10000LOCI_120D_01/samples \
               --output out \
               --npops 3 \
               --nloci 10000 \
               --pop_size 5000 \
               --seed 3168