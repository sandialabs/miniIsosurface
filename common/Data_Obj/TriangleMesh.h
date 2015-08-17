/*
 * TriangleMesh.h
 *
 *  Created on: Jul 7, 2015
 *      Author: sjmunn
 */

#ifndef DATA_OBJ_TRIANGLEMESH_H_
#define DATA_OBJ_TRIANGLEMESH_H_

#include"../includes.h"

#include"./Triplet.h"
#include"./MapReverse.h"

// Reporting Headers
#include"../Reporting/Log.h"

template<typename T>
class TriangleMesh {

	typedef Triplet<T> Point3d;
	typedef Triplet<T> Vector3d;
	typedef Triplet<unsigned> Triangle;
public:
	TriangleMesh();
	TriangleMesh(TriangleMesh&,MapReverse&);
	virtual ~TriangleMesh();

	// Data member interface: write
	void buildMesh(TriangleMesh&,MapReverse&);

	void addPoint(T *);
	void addPoint(T,T,T);
	void addPoint(Point3d);

	void addNormal(T *);
	void addNormal(T,T,T);
	void addNormal(Vector3d);

	void addTriangle(unsigned *);
	void addTriangle(unsigned,unsigned,unsigned);
	void addTriangle(Triangle);

	void setPoint(unsigned,const T *);
	void setNormal(unsigned,const T *);
	void resetTheMesh(void);

	// Combine 2 meshes
	TriangleMesh& operator+=(TriangleMesh&);

	// Data member interface: read
	const T* getPointPosition(unsigned) const;
	void getPointPosition(unsigned,Point3d &) const;
	const T* getNormalVector(unsigned) const;
	const unsigned* getTriangleIndices(unsigned) const;

	// Get info on the mesh
	unsigned numberOfVertices() const;
	unsigned numberOfTriangles() const;

	// This function is not a member b/c we keep I/O separate
	template <typename Anything> friend void saveTriangleMesh(const TriangleMesh<Anything> *, const char *);

private:
	// Prevent object copying
	TriangleMesh(const TriangleMesh&); // no implementation
	TriangleMesh& operator=(const TriangleMesh&); // no implementation

	std::vector<Point3d> points;
	std::vector<Vector3d> normals;
	std::vector<Triangle> indices;
};

#endif /* DATA_OBJ_TRIANGLEMESH_H_ */
