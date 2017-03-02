## openmpDupFree ##
When parallel versions are run, points along the boundaries may be accounted for multiple
times. This version is the same as openmp except it removes those duplicate points.
```
export OMP_NUM_THREADS=n
./openmp/openmp -i myImage.vtk -o outputMeshOpenMP.vtk -v 1.0 -g 1012
```

This version of the implementation uses openMP to speed up the algorithm. The image
is split up into uniform sections. Each openMP thread completes the algorithm over
a number of sections and then the result is combine to create one output mesh. Unlike
the `openmp` implementation, the output mesh of this implementation guarantees that
there are no duplicate points.

