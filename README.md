# Marching Cubes #

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
is done without any parallelization.


