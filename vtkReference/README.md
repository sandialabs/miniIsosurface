## VTK Reference Implementation ##

Example usage with an isovalue of 1.0:

```
./vtk/vtk -input_file myImage.vtk -output_file outputMeshSerial.vtk \
    -isoval 1.0 -algorithm mc
```

This is an implementation that uses VTK to run either marching cubes or
flying edges. Any parallelization that occurs comes from VTK. The VTK
library must be linked.

# Build Instructuions #

```
cmake /path/to/miniIsusurface/vtkReference -DVTK_DIR:PATH=/path/to/vtk/build
make
```


## License ##

miniIsosurface is distributed under the OSI-approved BSD 3-clause License.
See [LICENSE.txt](../LICENSE.txt) for details.

Copyright (c) 2017
National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under
the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
certain rights in this software.
