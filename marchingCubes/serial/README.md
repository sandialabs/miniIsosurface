## serial ##

Example usage with an isovalue of 1.0:

```
./serial/serial -input_file myImage.vtk -output_file outputMeshSerial.vtk \
    -isoval 1.0
```

This is the reference implementation. The algorithm is executed without any
parallelization.


## License ##

miniIsosurface is distributed under the OSI-approved BSD 3-clause License.
See [LICENSE.txt](../../LICENSE.txt) for details.

Copyright (c) 2017
National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under
the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
certain rights in this software.
