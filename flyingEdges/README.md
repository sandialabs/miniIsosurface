# miniIsosurface/flyingEdges #

The Flying Edges algorithm creates a polygonal mesh to approximate an isosurface
from a three-dimensional discrete scalar field. Similar to Marching Cubes, the
discrete scalar field is subdivided into uniform cubes. For each cube, the values
ateach corner are compared to the isovalue to determine a possible configuration of
triangles. Connected together, the triangles from all cubes will form the polygonal
mesh.

The difference between Flying Edges and Marching Cubes is that Flying Edges
traverses the grid row by row. By doing so, a small amount of additional work
reduces the total amoount of computation and the amount of memory required
for the output can be determined in advance.

# Build Instructions #
1. Clone the repository

```
git clone git@gitlab.sandia.gov:visMiniapps/miniIsosurface.git
```
2. Create a build directory. mantevo-marching-mubes does not require an out-of-tree build, but it is
cleaner.

```
mkdir miniIsosurface/feBuild
cd miniIsosurface/feBuild
```
3. Invoke CMake from your build directory, pointing to the miniIsosurface/flyingEdges source directory.
The flag BUILD\_CUDA is required to build the thrust implementation.

```
cmake /path/to/miniIsosurface/flyingEdges
```

or

```
cmake /path/to/miniIsosurface/flyingEdges -DBUILD_CUDA=ON
```

4. Invoke GNU make from the build directory.

```
make
```

After compiling, the following executables will be created
* `./serial/flyingEdgesSerial`
After compiling with BUILD\_CUDA, the following executables will also be created
* `./thrust/flyingEdgesThurst`

For all executables, the `input_file`, `output_file` and `isoval` flags must be set.
Upon running, each executable creates a yaml file describing performance and output
characteristics.

The contents of the yaml file is also printed to console.
To specify the yaml output file name, the flag is `yaml_output_file`.

Some executables have additional flags. To print out all flags for an executable,
use the `help` flag.
