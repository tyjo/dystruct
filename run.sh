#!/usr/bin/env bash
export OMP_NUM_THREADS=1

./bin/dystruct --input ./supp/data/BASELINE_100GEN_10000LOCI_120D_01/samples \
               --output out \
               --npops 3 \
               --nloci 10000 \
               --pop-size 5000 \
               --step-size-power -0.60 \
               --hold-out-fraction 0.2 \
               --seed 3168