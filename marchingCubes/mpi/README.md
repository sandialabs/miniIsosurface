## mpi ##
mpirun can be used to launch this executable. The -np flag tells mpi how many
processes to run.
```
mpirun -np numProcesses ./mpi/mpi -i myImage.vtk -o outputMeshMPI.vtk -v 1.0 -g 1012
```

This version of the implementation uses mpi to speed up the algorithm. The image
is split up into uniform sections. Each mpi process is in charge of reading it's
sections of the image from memory and running the algorihtm on those sections.
When reading sections from memory, the mpi process also reads in ghost cells.

There are two ways to output the mesh to file. The first is to output one mesh
to each process. In the example above, if `numProcesses=2` then two files would
be written, `outputMeshMPI.vtk.0` and `outputMeshMPI.vtk.1`. The other option
is to set the `one_mesh` flag to 1. In that case, once the algorithm has been
run and timed, the meshes from each process will be combined and one output
file will be written.


## License ##

miniIsosurface is distributed under the OSI-approved BSD 3-clause License.
See [LICENSE.txt](../../LICENSE.txt) for details.

Copyright (c) 2017
National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under
the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
certain rights in this software.
