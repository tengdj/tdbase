#!/bin/bash

for ((i=0; i<=9551; i++))
do
    ./tdbase join \
        --tile1 foo_n_nv50_nu200_vs100_r30.dt \
        --tile2 foo_v_nv50_nu200_vs100_r30.dt \
        -q nn \
        --knn 3 \
        -l 0 -l 20 -l 40 -l 60 -l 80 -l 100 \
        --cn 1 \
        --target $i
done