#!/bin/bash

#nn(nuclei-vessel)	gpu 100	gpu 50 100	aabb 100	aabb 50 100	multimbb 100	multimbb 50 100	

build/join --tile1 src/tmp_n_nv50_nu200_s10_vs400000.dt --tile2 src/tmp_v_nv50_nu200_s10_vs400000.dt -r 1000 --lod 100  --query nn --gpu >nn_nv_100_gpu 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs400000.dt --tile2 src/tmp_v_nv50_nu200_s10_vs400000.dt -r 1000 --lod 50 100  --query nn --gpu >nn_nv_50_100_gpu 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs400000.dt --tile2 src/tmp_v_nv50_nu200_s10_vs400000.dt -r 1000 --lod 100  --query nn --aabb >nn_nv_100_aabb 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs400000.dt --tile2 src/tmp_v_nv50_nu200_s10_vs400000.dt -r 1000 --lod 50 100  --query nn --aabb >nn_nv_50_100_aabb 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs100.dt --tile2 src/tmp_v_nv50_nu200_s10_vs100.dt -r 1000 --lod 100  --query nn >nn_nv_100_multibbb 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs100.dt --tile2 src/tmp_v_nv50_nu200_s10_vs100.dt -r 1000 --lod 50 100  --query nn >nn_nv_50_100_multibbb 2>&1

#nn(nuclei-nuclei)	gpu 100	gpu 50 100	aabb 100	aabb 50 100	
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs400000.dt -r 1000 --lod 100  --query nn --gpu >nn_nn_100_gpu 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs400000.dt -r 1000 --lod 50 100  --query nn --gpu >nn_nn_50_100_gpu 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs400000.dt -r 1000 --lod 100  --query nn --aabb >nn_nn_100_aabb 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs400000.dt -r 1000 --lod 50 100  --query nn --aabb >nn_nn_50_100_aabb 2>&1

#within 1000(nuclei-vessel)	gpu 100	gpu 50 100	aabb 100	aabb 50 100	multimbb 100	multimbb 50 100	multimbb 50 80 100	

build/join --tile1 src/tmp_n_nv50_nu200_s10_vs400000.dt --tile2 src/tmp_v_nv50_nu200_s10_vs400000.dt -r 1000 --query within --max_dist 1000 --lod 100 --gpu >within_nv_100_gpu 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs400000.dt --tile2 src/tmp_v_nv50_nu200_s10_vs400000.dt -r 1000 --query within --max_dist 1000 --lod 50 100 --gpu >within_nv_50_100_gpu 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs400000.dt --tile2 src/tmp_v_nv50_nu200_s10_vs400000.dt -r 1000 --query within --max_dist 1000 --lod 100 --aabb >within_nv_100_aabb 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs400000.dt --tile2 src/tmp_v_nv50_nu200_s10_vs400000.dt -r 1000 --query within --max_dist 1000 --lod 50 100  --aabb >within_nv_50_100_aabb 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs100.dt --tile2 src/tmp_v_nv50_nu200_s10_vs100.dt -r 1000 --query within --max_dist 1000 --lod 100 >within_nv_100_multimbb 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs100.dt --tile2 src/tmp_v_nv50_nu200_s10_vs100.dt -r 1000 --query within --max_dist 1000 --lod 50 100 >within_nv_50_100_multimbb 2>&1


#within 1000(nuclei-nuclei)	gpu 100	gpu 50 100	aabb 100	aabb 50 100	
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs400000.dt -r 1000 --lod 100  --query within --gpu >within_nn_100_gpu 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs400000.dt -r 1000 --lod 50 100  --query within --gpu >within_nn_50_100_gpu 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs400000.dt -r 1000 --lod 100  --query within --aabb >within_nn_100_aabb 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs400000.dt -r 1000 --lod 50 100  --query within --aabb >within_nn_50_100_aabb 2>&1
