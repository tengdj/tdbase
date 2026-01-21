# TDBase
we propose a novel multi-level progressive 3D compression method, which is designed to satisfy two progressive query conditions, thus supports returning correct results even with compressed data whenever possible. Based on this, we further propose novel progressive refinement which starts from lightweight low-resolution representation for possible early return of results and progresses to higher resolution representation, which can be largely avoided. A system TDBase is implemented with progressive refinement to support highly efficient and accurate spatial queries for complex three-dimensional objects.

## install and setup
install the prerequisite libraries: cmake CGAL GMP Boost(program_options) ZLIB OpenMP
if the machine equiped with NVIDIA GPU, setup the environment before compiling
```console
export USE_GPU=TRUE_

```
then compile the code
```console
cd src/
mkdir build
cd build/
cmake ../
make simulator
make tdbase

```

## generate synthetic dataset
We developped a tool to generate synthetic datasets 
```console
Simulator:
  -h, --help                     produce help message
  --hausdorff                    enable Hausdorff distance calculation
  -m, --multi_lods               the input are polyhedrons in multiple files
  -i, --allow_intersection       allow the nuclei to intersect with other nuclei or vessel
  -n, --nuclei arg (=nuclei.pt)  path to the nuclei prototype file
  -v, --vessel arg (=vessel.pt)  path to the vessel prototype file
  -o, --output arg (=default)    prefix of the output files
  -t, --threads arg (=8)         number of threads
  --nv arg (=50)                 number of vessels
  --nu arg (=200)                number of nucleis per vessel
  --vs arg (=100)                number of vertices in each voxel
  --verbose arg (=0)             verbose level
  -r, --sample_rate arg (=30)    sampling rate for Hausdorff distance calculation

```
You can generate datasets with varying sizes and configurations. The below command generated two .dt files which contain 50 complex vessels and about 10000 (50*200) nuclei. queries can be conducted over those two dt files foo_n_nv50_nu200_vs100_r30.dt and foo_v_nv50_nu200_vs100_r30.dt. 
```console
./simulator -n ../../data/nuclei.pt -v ../../data/vessel.pt -o foo --hausdorff --nv 50 --nu 200
```
Note that for now we disabled the facet-association-based optimization for hausdorff distance calculation but depend solely on AABB-tree, thus it takes a while to complete the generation with hausdorff distances calculated. In addition, you can use the pre-generated synthetic dt files in the data folder. 

The generated files are in compressed format which need to be decoded during querying. Converting .dt files into decoded format can achieve better query performance but with significant higher storage cost. We are working on improving the decoding efficiency with partial decoding and GPU acceleration.

```console
./tdbase convert foo_n_nv50_nu200_vs100_r30.dt decoded_n_nv50_nu200_vs100_r30.dt
./tdbase convert foo_v_nv50_nu200_vs100_r30.dt decoded_v_nv50_nu200_vs100_r30.dt

```

## run test
Three queries are implemented: KNN, within-distance, and intersection. For now, the intersection query is still buggy and we are working on it, currently it can be tested using the within-distance query by setting the distance threshold to 0. 
```console
-h, --help                          produce help message
  -t, --threads arg (=8)              number of threads
  --cn arg (=1)                       number of threads for geometric computation for each tile
  -g, --gpu                           compute with GPU
  -v, --verbose arg (=0)              verbose level
  -p, --print_result                  print result to standard out
  --aabb                              calculate distance with aabb
  -c, --counter_clock                 is the faces recorded clock-wise or counterclock-wise
  --disable_byte_encoding             using the raw hausdorff instead of the byte encoded ones
  -l, --lods arg                      the lods that needs be processed
  --tile1 arg                         path to tile 1
  --tile2 arg                         path to tile 2
  --max_objects1 arg (=9223372036854775807)
                                      max number of objects in tile 1
  --max_objects2 arg (=9223372036854775807)
                                      max number of objects in tile 2
  -q, --query arg                     query types: intersect|nn|within
  -k, --knn arg (=1)                  the K value for NN query
  -w, --within_dist arg (=1000)       the maximum distance for within query
```

The below command conduct a 3NN join. Specifying "-g" option makes the geometry computation conducted with GPU. In cases GPU is not supported, you can specify how many threads can be used for geometric computations with option "--cn 10" (using 10 threads), the default value is 1. 
```console
./tdbase join --tile1 decoded_n_nv50_nu200_vs100_r30.dt --tile2 decoded_v_nv50_nu200_vs100_r30.dt -q nn --knn 3 -l 20 -l 40 -l 60 -l 80 -l 100 --cn 10 -p > result.txt

```

The "-p" option enabling the printing of the final results to the standard out, and stored in result.txt. You can also run test without progressive querying by querying directly on the highest LOD with the following command. 
```console
./tdbase join --tile1 decoded_n_nv50_nu200_vs100_r30.dt --tile2 decoded_v_nv50_nu200_vs100_r30.dt -q nn --knn 3 -l 100 --cn 10 -p > result100.txt

```

Changing "-l 100" to other LODs like "-l 20" can query specific LOD and obtain approximate result. We developped a tool to compare the result with the ground truth (the result while querying LOD 100). 

```console
./tdbase evaluate result100.txt result.txt

```

The following command conducts a 3NN query on a single dataset as a self-join.  
```console
./tdbase join --tile1 decoded_n_nv50_nu200_vs100_r30.dt -q nn --knn 3 -l 20 -l 40 -l 60 -l 80 -l 100 --cn 10 -p > result.txt

```

The following command conducts a within distance query, and the distance threshold is set to 100. 
```console
./tdbase join --tile1 foo_n_nv50_nu200_vs100_r30.dt --tile2 foo_v_nv50_nu200_vs100_r30.dt -q within -w 100 -l 20 -l 40 -l 60 -l 80 -l 100 -p >result.txt

```

You can also conduct a within-distance query without progressive query on the highest LOD, and evaluate the accuracy as discussed above. 




