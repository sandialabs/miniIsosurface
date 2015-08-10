/*
 * TriangleMesh.cpp
 *
 *  Created on: Jul 7, 2015
 *      Author: sjmunn
 */

#include "TriangleMesh.h"

// Constructor/destructors

template<typename T>
TriangleMesh<T>::TriangleMesh() {
	// Not needed
}

template<typename T>
TriangleMesh<T>::~TriangleMesh() {
	// Not needed
}

// Data member interface: write

template<typename T>
void TriangleMesh<T>::addPoint(T * coordinates) {
	Point3d pt(coordinates);
	points.push_back(pt);
}

template<typename T>
void TriangleMesh<T>::addPoint(T a, T b, T c) {
	Point3d pt(a,b,c);
	points.push_back(pt);
}

template<typename T>
void TriangleMesh<T>::addPoint(Point3d pt) {
	points.push_back(pt);
}

template<typename T>
void TriangleMesh<T>::addNormal(T * normal) {
	Vector3d norm(normal);
	normals.push_back(norm);
}

template<typename T>
void TriangleMesh<T>::addNormal(T a, T b, T c) {
	Vector3d norm(a,b,c);
	normals.push_back(norm);
}

template<typename T>
void TriangleMesh<T>::addNormal(Vector3d norm) {
	normals.push_back(norm);
}


template<typename T>
void TriangleMesh<T>::addTriangle(unsigned * triIndices) {
	Triangle TriIndex(triIndices);
	indices.push_back(TriIndex);
}

template<typename T>
void TriangleMesh<T>::addTriangle(unsigned a, unsigned b, unsigned c) {
	Triangle TriIndex(a,b,c);
	indices.push_back(TriIndex);
}

template<typename T>
void TriangleMesh<T>::addTriangle(Triangle TriIndex) {
	indices.push_back(TriIndex);
}

template<typename T>
TriangleMesh<T>& TriangleMesh<T>::operator+=(TriangleMesh<T>& newMesh) {
	unsigned nPointsStart = this->numberOfVertices();
	unsigned nTrianglesNew = newMesh.numberOfTriangles();

	this->points.insert(this->points.end(), newMesh.points.begin(),newMesh.points.end());
	this->normals.insert(this->normals.end(), newMesh.normals.begin(),newMesh.normals.end());

	for (unsigned iTriangle=0;iTriangle<nTrianglesNew;iTriangle++) {
		newMesh.indices[iTriangle]+=nPointsStart;
	}
	this->indices.insert(this->indices.end(), newMesh.indices.begin(), newMesh.indices.end());

	return *this;
}

// Data member interface: read

template<typename T>
const T* TriangleMesh<T>::getPointPosition(unsigned idx) const {
	return points[idx].getCoordinates();
}

template<typename T>
void TriangleMesh<T>::getPointPosition(unsigned idx, Point3d& returnPoint) const {
	returnPoint=points[idx];
}

template<typename T>
const T* TriangleMesh<T>::getNormalVector(unsigned idx) const {
	return normals[idx].getCoordinates();
}

template<typename T>
const unsigned* TriangleMesh<T>::getTriangleIndices(unsigned idx) const {
	return indices[idx].getCoordinates();
}

// Get info on the mesh

template<typename T>
unsigned TriangleMesh<T>::numberOfVertices() const {
	return points.size();
}

template<typename T>
unsigned TriangleMesh<T>::numberOfTriangles() const {
	return indices.size();
}

// Must instantiate class for separate compilation
template class TriangleMesh<float_t>;
