/*
 * SaveTriangleMesh.h
 *
 *  Created on: Jul 10, 2015
 *      Author: sjmunn
 */

#ifndef SAVETRIANGLEMESH_H_
#define SAVETRIANGLEMESH_H_

#include"../includes.h"

#include "ConvertBuffer.h"

// Data Objects
#include"../types.h"

// This Function is a friend of the Triangle Mesh class and has access to private members of that class
template<typename T>
void saveTriangleMesh(const TriangleMesh<T> *mesh, const char *vtkFileName) {
	std::ofstream stream;
	stream.open(vtkFileName);
	int spacialDimensions = 3;

	unsigned nverts = mesh->numberOfVertices();
	unsigned ntriangles = mesh->numberOfTriangles();
	TypeInfo ti = createTemplateTypeInfo<T>();

	size_t bufsize = 0;
	std::vector<char> wbuff;

	stream << "# vtk DataFile Version 3.0" << std::endl;
	stream << "Isosurface Mesh" << std::endl;
	stream << "BINARY" << std::endl;
	stream << "DATASET POLYDATA" << std::endl;

	bufsize = mesh->points.size() * spacialDimensions * sizeof(T);
	wbuff.resize(bufsize);

	// Writing points data
	stream << "POINTS " << nverts << " " << ti.name() << std::endl;

	T *bufPointer = reinterpret_cast<T*>(&wbuff[0]);

	typedef Triplet<T> Point3d;
	for (typename std::vector<Point3d>::const_iterator it=mesh->points.begin(); it != mesh->points.end(); ++it) {
		for (int iCoordinate = 0; iCoordinate < 3; ++iCoordinate) {
			*bufPointer = (*it)[iCoordinate];
			flipEndianness(*bufPointer++);
		}
	}

	stream.write(&wbuff[0], wbuff.size());
	stream << std::endl;

	// Writing triangle indices

	bufsize = ntriangles * 4 * sizeof(unsigned);
	wbuff.resize(bufsize);
	unsigned *ind = reinterpret_cast<unsigned*>(&wbuff[0]);

	typedef Triplet<unsigned> Triangle;
	for (typename std::vector<Triangle>::const_iterator it=mesh->indices.begin(); it != mesh->indices.end(); ++it) {
		*ind = 3;
		flipEndianness(*ind++);
		for (int j = 0; j < 3; ++j) {
			*ind = (*it)[j];
			flipEndianness(*ind++);
		}
	}

	stream << "POLYGONS " << ntriangles << " " << ntriangles * 4 << std::endl;
	stream.write(&wbuff[0], wbuff.size());
	stream << std::endl;

	bufsize = mesh->normals.size() * spacialDimensions * sizeof(T);
	wbuff.resize(bufsize);
	bufPointer = reinterpret_cast<T*>(&wbuff[0]);

	for (typename std::vector<Point3d>::const_iterator it=mesh->normals.begin(); it != mesh->normals.end(); ++it) {
		for (int iCoordinate = 0; iCoordinate < 3; ++iCoordinate) {
			*bufPointer = (*it)[iCoordinate];
			flipEndianness(*bufPointer++);
		}
	}

	stream << "POINT_DATA " << nverts << std::endl;
	stream << "NORMALS Normals " << ti.name() << std::endl;
	stream.write(&wbuff[0], wbuff.size());
	stream << std::endl;

	stream.close();
}

#endif /* SAVETRIANGLEMESH_H_ */
