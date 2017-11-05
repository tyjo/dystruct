#!/usr/bin/env bash

for folder in $(find ../data ! -name "simulate_data.sh" -type d);
do
    ./format.py $folder
done