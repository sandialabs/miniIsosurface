# SameContentsCheck #

There are multiple ways that an output mesh can be represented. To check
that two output files from mantevo-marching-cubes have equivalent meshes,
use SameContentsCheck:

```
./tests/SameContentsCheck outputMeshSerial.vtk outputMeshOpenMP.vtk
```


## License ##

miniIsosurface is distributed under the OSI-approved BSD 3-clause License.
See [LICENSE.txt](../../LICENSE.txt) for details.

Copyright (c) 2017
National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under
the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
certain rights in this software.
