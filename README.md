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
* The `std::to_string` function

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
The two flags, BUILD\_OPENMP and BUILD\_MPI, must be turned on to build all implementations.

```
cmake /path/to/mantevo-marching-cubes -DBUILD_OPENMP=On -DBUILD_MPI=On
```

4. Invoke GNU make from the build directory.

```
make
```

After compiling, the following executables will be created:
* `./serial/serial`
* `./openmp/openmp`
* `./openmpDupFree/openmpDupFree`
* `./mpi/mpi`
* `./openmpAndMpi/openmpAndMpi`
* `./tests/SameContentsCheck`

For all executables, the `input_file`, `output_file` and `isoval` flags must be set.
Upon running, each executable creates a yaml file describing performance and output
characteristics.

The contents of the yaml file is also printed to console.
To specify the yaml output file name, the flag is `yaml_output_file`.

Some executables have additional flags. To print out all flags for an executable,
use the `help` flag.

## Kokkos Build Instructions ##

The Kokkos implementation does not use CMake. See the Kokkos README for the Marching Cubes
with Kokkos specific build instrucions.
