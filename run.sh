#!/usr/bin/env bash
export OMP_NUM_THREADS=2

./bin/dystruct --input ./supp/example_data/samples.geno \
               --generation-times ./supp/example_data/sample_times \
               --output out \
               --npops 3 \
               --nloci 10000 \
               --pop-size 5000 \
               --seed 1145 \
               --hold-out-fraction 0.1 \
               --hold-out-seed 55307