# openmpAndMpi #

mpirun can be used to launch this executable. The -np flag tells mpi how
many processes to run.

```
mpirun -np numProcesses ./openmpAndMpi/openmpAndMpi -i myImage.vtk \
    -o outputMeshOpenmpAndMpi -v 1.0 -g 1012
```

This implementation is a combination of the `openmpDupFree` and `mpi`
implementations. It uses both mpi and openmp to speed up the algorithm.


## License ##

miniIsosurface is distributed under the OSI-approved BSD 3-clause License.
See [LICENSE.txt](../../LICENSE.txt) for details.

Copyright (c) 2017
National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under
the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
certain rights in this software.
