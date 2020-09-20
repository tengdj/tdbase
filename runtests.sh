#!/bin/bash

#nn(nuclei-vessel)	gpu 100	gpu 50 100	aabb 100	aabb 50 100	multimbb 100	multimbb 50 100	

build/join --tile1 src/tmp_n_nv50_nu200_s10_vs100.dt --tile2 src/tmp_v_nv50_nu200_s10_vs100.dt -r 1000 --query nn --lod 100 --gpu >nn_nv_100_gpu.out 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs100.dt --tile2 src/tmp_v_nv50_nu200_s10_vs100.dt -r 1000 --query nn --lod 50 100 --gpu >nn_nv_50_100_gpu.out 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs100.dt --tile2 src/tmp_v_nv50_nu200_s10_vs100.dt -r 95 --query nn --lod 100 --aabb >nn_nv_100_aabb.out 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs100.dt --tile2 src/tmp_v_nv50_nu200_s10_vs100.dt -r 1000 --query nn --lod 50 100 --aabb >nn_nv_50_100_aabb.out 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs100.dt --tile2 src/tmp_v_nv50_nu200_s10_vs100.dt -r 1000 --query nn --lod 100 --multiple_mbb >nn_nv_100_multibbb.out 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs100.dt --tile2 src/tmp_v_nv50_nu200_s10_vs100.dt -r 1000 --query nn --lod 50 100 --multiple_mbb >nn_nv_50_100_multibbb.out 2>&1

#nn(nuclei-nuclei)	gpu 100	gpu 50 100	aabb 100	aabb 50 100	
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs100.dt -r 1000 --lod 100  --query nn --gpu >nn_nn_100_gpu.out 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs100.dt -r 1000 --lod 50 100  --query nn --gpu >nn_nn_50_100_gpu.out 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs100.dt -r 1000 --lod 100  --query nn --aabb >nn_nn_100_aabb.out 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs100.dt -r 1000 --lod 50 100  --query nn --aabb >nn_nn_50_100_aabb.out 2>&1

#within 1000(nuclei-vessel)	gpu 100	gpu 50 100	aabb 100	aabb 50 100	multimbb 100	multimbb 50 100	multimbb 50 80 100	

build/join --tile1 src/tmp_n_nv50_nu200_s10_vs100.dt --tile2 src/tmp_v_nv50_nu200_s10_vs100.dt -r 1000 --query within --max_dist 1000 --lod 100 --gpu >within_nv_100_gpu.out 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs100.dt --tile2 src/tmp_v_nv50_nu200_s10_vs100.dt -r 1000 --query within --max_dist 1000 --lod 50 100 --gpu >within_nv_50_100_gpu.out 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs100.dt --tile2 src/tmp_v_nv50_nu200_s10_vs100.dt -r 1000 --query within --max_dist 1000 --lod 100 --aabb >within_nv_100_aabb.out 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs100.dt --tile2 src/tmp_v_nv50_nu200_s10_vs100.dt -r 1000 --query within --max_dist 1000 --lod 50 100  --aabb >within_nv_50_100_aabb.out 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs100.dt --tile2 src/tmp_v_nv50_nu200_s10_vs100.dt -r 1000 --query within --max_dist 1000 --lod 100 --multiple_mbb >within_nv_100_multimbb.out 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs100.dt --tile2 src/tmp_v_nv50_nu200_s10_vs100.dt -r 1000 --query within --max_dist 1000 --lod 50 100 --multiple_mbb >within_nv_50_100_multimbb.out 2>&1


#within 1000(nuclei-nuclei)	gpu 100	gpu 50 100	aabb 100	aabb 50 100	
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs100.dt -r 1000 --query within --lod 100 --gpu >within_nn_100_gpu.out 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs100.dt -r 1000 --query within --lod 50 100  --gpu >within_nn_50_100_gpu.out 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs100.dt -r 1000 --query within --lod 100 --aabb >within_nn_100_aabb.out 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs100.dt -r 1000 --query within --lod 50 100 --aabb >within_nn_50_100_aabb.out 2>&1

#intersect	100	50 100
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs100.dt --tile2 src/tmp_n2_nv50_nu200_s10_vs100.dt -r 1000 --query intersect --lod 100 >intersect_nn_100.out 2>&1
build/join --tile1 src/tmp_n_nv50_nu200_s10_vs100.dt --tile2 src/tmp_n2_nv50_nu200_s10_vs100.dt -r 1000 --query intersect --lod 50 100 >intersect_nn_50_100.out 2>&1