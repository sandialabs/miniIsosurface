/* Generate synthetic data
   for marching cubes */

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>

inline void flipEndianness(float &val) {
	int mid = sizeof(float) / 2;
	char *data = reinterpret_cast<char*>(&val);
	for (int i = 0; i < mid; ++i) {
		std::swap(data[i], data[sizeof(float) - i - 1]);
	}
}

int main(void) {
  std::ofstream stream("dataGen.vtk",std::ofstream::binary);
  
  const unsigned dim = 1024;
  const unsigned long nPoints = dim*dim*dim;

  const float mins[3] = { -1, -1, -1 };
  const float maxs[3] = { 1, 1, 1 };
  
  float data[dim]; // Save full row of points at a time
  int writeSize = dim * sizeof(float);

  /* Write the file header */
  stream << "# vtk DataFile Version 3.0" << std::endl;
  stream << "vtk output" << std::endl;
  stream << "BINARY" << std::endl;
  stream << "DATASET STRUCTURED_POINTS" << std::endl;
  stream << "DIMENSIONS " << dim << " " << dim << " " << dim  << std::endl;
  stream << "SPACING " << "1 1 1" << std::endl;
  stream << "ORIGIN 0 0 0" << std::endl;
  stream << "POINT_DATA " << nPoints << std::endl;
  stream << "SCALARS ImageFile float" << std::endl;
  stream << "LOOKUP_TABLE default" << std::endl;

  // Port of Kewei's data generation code from his isosurface implemention
  // See isosurface.cpp:TangleField
  for (int z = 0; z < dim; ++z){
    for (int y = 0; y < dim; ++y){
      for (int x = 0; x < dim; ++x){
        const float xx = 3.0 * (mins[0] + (maxs[0] - mins[0]) * (1.0 * x / dim));
        const float yy = 3.0 * (mins[1] + (maxs[1] - mins[1]) * (1.0 * y / dim));
        const float zz = 3.0 * (mins[2] + (maxs[2] - mins[2]) * (1.0 * z / dim));
        const float v = (xx*xx*xx*xx - 5.0f*xx*xx + yy*yy*yy*yy - 5.0f*yy*yy
			 + zz*zz*zz*zz - 5.0f*zz*zz + 11.8f) * 0.2f + 0.5f;
	data[x]=v;
	flipEndianness(data[x]);
      }
      
      stream.write((char *)data,writeSize);
    }
  }
}
