# mantevo-marching-cubes #

Marching Cubes is an algorithm that creates a polygonal mesh to approximate an isosurface
from a three-dimensional discrete scalar field. The discrete scalar field is subdivided into
uniform cubes. For each cube, the values at each corner are compared to the isovalue
to determine a possible configuration of triangles. Connected together, the triangles from all
cubes will form the polygonal mesh.

Each vertex of a triangle is located on an edge of the cube. The location is approximated
by linearly interpoloating from the two endpoints of that edge. Similarly, the normal value of the
surface located at each vertex is approximated by linerly interpolating from the normals located at the
two endpoints of the edge where the normals are equal to the gradient, approximated by the difference
method.

This implementation of Marching Cubes has several implementations. The reference implementation
is done without any parallelization. The openmp implementation uses openMP.

# C++11 usage #

mantevo-marching-cubes makes use of C++11 features. Namely
* `auto` type-specifier
* range based for-loops
* The `std::vector<T>::data` member function
* The `std::vector<T>::emplace_back` member function
* The `std::array<T, N>` class

# Build Instructions #
1. Clone the repository

```
git clone https://gitlab.sandia.gov/kmorel-src/mantevo-marching-cubes.git
```
2. Create a build directory. mantevo-marching-mubes does not require an out-of-tree build, but it is
cleaner.

```
mkdir my_build_dir
cd my_build_dir
```
3. Invoke CMake from your build directory, pointing to the mantevo-marching-cubes source directory.

```
cmake /path/to/mantevo-marching-cubes
```
4. Invoke GNU make from the build directory.

```
make
```

After compiling, the following executables will be created:
* `./serial/serial`
* `./openmp/openmp`
* `./openmpDupFree/openmpDupFree`
* `./tests/SameContentsCheck`

## serial ##
Example usage with an isovalue of 1.0:
```
./serial/serial myImage.vtk outputMeshSerial.vtk 1.0
```

## openmp ##
The openmp executable contains an extra argument from the reference implementation, the
grianDim size. The grainDim size controls the granularity of the parallel work.
Example usage with an isovalue of 1.0 and grainDim size of 1012:
```
export OMP_NUM_THREADS=n
./openmp/openmp myImage.vtk outputMeshOpenMP.vtk 1.0 1012
```

## openmpDupFree ##
When parallel versions are run, points along the boundaries may be accounted for multiple
times. This version is the same as openmp except it removes those duplicate points.
```
export OMP_NUM_THREADS=n
./openmpDupFree/openmpDupFree myImage.vtk outputMeshOpenMP_NoSame.vtk 1.0 1012
```

## SameContentsCheck ##
There are multiple ways that an output mesh can be represented. To check that two output files
from mantevo-marching-cubes have equivalent meshes, use SameContentsCheck:
```
./tests/SameContentsCheck outputMeshSerial.vtk outputMeshOpenMP.vtk
```

