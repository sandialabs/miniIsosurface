## openmp ##

The openmp executable contains an extra argument from the reference
implementation, the grianDim size. The grain_dim size controls the
granularity of the parallel work. Example usage with an isovalue of 1.0 and
grain_dim size of 1012:

```
export OMP_NUM_THREADS=n
./openmp/openmp -input_file myImage.vtk -output_file outputMeshOpenMP.vtk \
    -isoval 1.0 -grain_dim 1012
```

This version of the implementation uses openMP to speed up the algorithm.
The image is split up into uniform sections. Each openMP thread completes
the algorithm over a number of sections and then the result is combine to
create one output mesh.

If two openMP threads have triangles that contain a point on the same edge,
then that point will appear twice in the output mesh. The `openmpDupFree`
implementation removes duplicate points.


## License ##

miniIsosurface is distributed under the OSI-approved BSD 3-clause License.
See [LICENSE.txt](../../LICENSE.txt) for details.

Copyright (c) 2017
National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under
the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
certain rights in this software.
