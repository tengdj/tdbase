# TDBase
we propose a novel multi-level progressive 3D compression method, which is designed to satisfy two progressive query conditions, thus supports returning correct results even with compressed data whenever possible. Based on this, we further propose novel progressive refinement which starts from lightweight low-resolution representation for possible early return of results and progresses to higher resolution representation, which can be largely avoided. A system TDBase is implemented with progressive refinement to support highly efficient and accurate spatial queries for complex three-dimensional objects.
## install and setup
install the prerequisite libraries: cmake CGAL GMP Boost(program_options) ZLIB OpenMP

```console
cd src/
mkdir build
cd build/
cmake ../
make simulator
make tdbase

```

## generate synthetic dataset
generate two .dt files, one for nuclei (foo_n_nv1000_nu200_vs100_r30_cm1.dt) and one for vessel (foo_v_nv1000_nu200_vs100_r30_cm1.dt). --nv specifies the number of vessels in the dataset while --nu specifies the number of nuclei around each vessel. Thus the generated dataset contains 1000 vessels and 200,000 nuclei
```console
./simulator -u ../../data/nuclei.pt -v ../../data/vessel.pt -o foo --nv 1000 --nu 200

```

## run test
conduct a 3NN join. -g specifies that the geometry computation are conducted with GPU. 
```console
./tdbase join --tile1 foo_n_nv1000_nu200_vs100_r30_cm1.dt --tile2 foo_v_nv1000_nu200_vs100_r30_cm1.dt -q nn --knn 3 -g --lod 20 40 60 80 100

```


conduct a 3NN join with the progressive refinement disabled by checking only the highest LOD polyhedrons.
```console
./tdbase join --tile1 foo_n_nv1000_nu200_vs100_r30_cm1.dt --tile2 foo_v_nv1000_nu200_vs100_r30_cm1.dt -q nn --knn 3 -g --lod 100

```

conduct a 3NN join with the geometric computation conducted using 24 threads (--cn 24)
```console
./tdbase join --tile1 foo_n_nv1000_nu200_vs100_r30_cm1.dt --tile2 foo_v_nv1000_nu200_vs100_r30_cm1.dt -q nn --knn 3 --cn 24 --lod 20 40 60 80 100

```
conduct a within distance join which conducts a within 50 distance join.
```console
./tdbase join --tile1 foo_n_nv1000_nu200_vs100_r30_cm1.dt --tile2 foo_v_nv1000_nu200_vs100_r30_cm1.dt -q within --within_dist 50 -g --lod 20 40 60 80 100 

```
conduct a 3NN self join.
```console
./tdbase join --tile1 foo_n_nv1000_nu200_vs100_r30_cm1.dt -q nn --knn 3 -g --lod 20 40 60 80 100

```


