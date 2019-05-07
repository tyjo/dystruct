#!/usr/bin/env bash

for i in $(seq 2 3); do
    python plot_Q.py ./plot_demo/Q${i} \
                     ./plot_demo/samplelabels \
                     ./plot_demo/poporder \
                     --match-Q ./plot_demo \
                     --spacing 8 \
                     --width 8 \
                     --height 3 \
                     --fontsize 8
    mv dystruct_k${i}.pdf plot_demo/dystruct_k${i}.pdf
done