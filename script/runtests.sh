#!/bin/bash

#intersect

	# nuclei
		#brute
#../build/join -q intersect --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt --tile2 teng_n2_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 100
#../build/join -q intersect --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt --tile2 teng_n2_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 20 100

# nearest neighbor
	# nuclei
		# brute
../build/join -q nn --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 100
../build/join -q nn --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 40 60 100
		# partition
../build/join -q nn --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 100 -m
../build/join -q nn --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 40 60 100 -m
		# AABB
../build/join -q nn --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 100 --aabb
../build/join -q nn --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 40 60 100 --aabb
		# GPU
../build/join -q nn --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 100 -g
../build/join -q nn --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 40 60 100 -g

	# vessel
		#brute
../build/join -q nn --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt --tile2 teng_v_nv50_nu200_s10_vs100_r30.dt -n 50 -r 100 --lod 100
../build/join -q nn --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt --tile2 teng_v_nv50_nu200_s10_vs100_r30.dt -n 50 -r 100 --lod 40 100
		#partition
../build/join -q nn --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt --tile2 teng_v_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 100 -m
../build/join -q nn --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt --tile2 teng_v_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 40 60 100 -m
		#AABB
../build/join -q nn --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt --tile2 teng_v_nv50_nu200_s10_vs100_r30.dt -n 50 -r 100 --lod 100 --aabb
../build/join -q nn --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt --tile2 teng_v_nv50_nu200_s10_vs100_r30.dt -n 50 -r 100 --lod 40 100 --aabb
		#GPU
../build/join -q nn --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt --tile2 teng_v_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 100 -g
../build/join -q nn --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt --tile2 teng_v_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 40 100 -g
		
# within
	# nuclei
		# brute
../build/join -q within --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 100
../build/join -q within --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 40 60 100
		# partition
../build/join -q within --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 100 -m
../build/join -q within --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 40 60 100 -m
		# AABB
../build/join -q within --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 100 --aabb
../build/join -q within --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 40 60 100 --aabb
		# GPU
../build/join -q within --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 100 -g
../build/join -q within --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 40 60 100 -g

	# vessel
		#brute
../build/join -q within --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt --tile2 teng_v_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 100
../build/join -q within --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt --tile2 teng_v_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 50 60 70 100
		#partition
../build/join -q within --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt --tile2 teng_v_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 100 -m
../build/join -q within --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt --tile2 teng_v_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 50 60 70 100 -m
		#AABB
../build/join -q within --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt --tile2 teng_v_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 100 --aabb
../build/join -q within --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt --tile2 teng_v_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 50 60 70 100 --aabb
		#GPU
../build/join -q within --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt --tile2 teng_v_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 100 -g
../build/join -q within --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt --tile2 teng_v_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 50 60 70 100 -g
