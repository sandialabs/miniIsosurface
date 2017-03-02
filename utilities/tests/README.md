# SameContentsCheck #
There are multiple ways that an output mesh can be represented. To check that two output files
from mantevo-marching-cubes have equivalent meshes, use SameContentsCheck:
```
./tests/SameContentsCheck outputMeshSerial.vtk outputMeshOpenMP.vtk
```
