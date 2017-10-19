# Marching Cubes with Kokkos Build Instructions #

1. Go to the kokkos folder within mantevo-marching-cubes
```
cd mantevo-marching-cubes/build/dir/
cd kokkos/
```
2. Invoke make with path to the Kokkos directory
```
make -e KOKKOS_PATH=/Path/To/kokkos/
```
After compiling, the following executable will be created:
* `./kokkos.host`

# kokkos #

This implementation if similar to the openmpDupFree implementation with
the exception that all memory allocation for the parallel portion of the
code is done ahead of time. In fact, if not enough memory is available,
an exception will be thrown. The number of points and triangles allocated
can be controlled by the command line flags `points_allocate` and
`triangles_allocate`. Left unset, the algorithm will make a guess as to
how much to allocate.

```
./kokkos -i myImages.vtk -o outputMeshKokkos.vtk -v 1.0 -g 1012 \
    -points_allocate 2000000 -triangles_allocate 4000000
```


## License ##

miniIsosurface is distributed under the OSI-approved BSD 3-clause License.
See [LICENSE.txt](../../LICENSE.txt) for details.

Copyright (c) 2017
National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under
the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
certain rights in this software.
