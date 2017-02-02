# openmpAndMpi #
mpirun can be used to launch this executable. The -np flag tells mpi how many
processes to run.
```
mpirun -np numProcesses ./openmpAndMpi/openmpAndMpi -i myImage.vtk \
  -o outputMeshOpenmpAndMpi -v 1.0 -g 1012
```

This implementation is a combination of the `openmpDupFree` and `mpi` implementations.
It uses both mpi and openmp to speed up the algorithm.
