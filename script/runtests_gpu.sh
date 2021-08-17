#!/bin/bash
# nearest neighbor
	# nuclei
		# GPU
../build/join -q nn --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 100 -g
../build/join -q nn --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 40 60 100 -g

	# vessel
		#GPU
../build/join -q nn --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt --tile2 teng_v_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 100 -g
../build/join -q nn --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt --tile2 teng_v_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 40 100 -g
		
# within
	# nuclei
		# GPU
../build/join -q within --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 100 -g
../build/join -q within --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 40 60 100 -g

	# vessel
		#GPU
../build/join -q within --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt --tile2 teng_v_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 100 -g
../build/join -q within --tile1 teng_n_nv50_nu200_s10_vs100_r30.dt --tile2 teng_v_nv50_nu200_s10_vs100_r30.dt -n 50 -r 1000 --lod 50 60 70 100 -g
