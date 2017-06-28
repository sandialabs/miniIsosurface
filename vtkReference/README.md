## VTK Reference Implementation ##
Example usage with an isovalue of 1.0:
```
./vtk/vtk -input_file myImage.vtk -output_file outputMeshSerial.vtk -isoval 1.0 -algorithm mc
```

This is an implementation that uses VTK to run either marching cubes or flying
edges. Any parallelization that occurs comes from VTK. The VTK library must be
linked.

# Build Instructuions #

```
cmake /path/to/miniIsusurface/vtkReference -DVTK_DIR:PATH=/path/to/vtk/build
make
```
