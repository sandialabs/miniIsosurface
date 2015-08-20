/*
 * buildMesh.h
 *
 *  Created on: Aug 19, 2015
 *      Author: sjmunn
 */

#ifndef ALGORITHM_BUILDMESH_H_
#define ALGORITHM_BUILDMESH_H_

// External includes
#include"../../common/includes.h"

// Reporting Headers
#include"../../common/Reporting/YAML_Element.hpp"
#include"../../common/Reporting/YAML_Doc.hpp"
#include"../../common/Reporting/Log.h"

// Strategy base class
#include"../../common/Data_Obj/TriangleMesh.h"

#include"./DuplicateRemover.h"

template<typename T>
void buildMesh(TriangleMesh<T>& newMesh, TriangleMesh<T>& ogMesh,DuplicateRemover& newPtMap) {
	/*
	 * Iterates through the entire dataArray in newPtMap and creates
	 * a new mesh from the new point indices in newPointIdx
	 */
	unsigned nOgPoints=newPtMap.getSize();
	CLOG(logDEBUG) << "nogpoints: " << nOgPoints;
	// As many new points as last index in the dataArray
	unsigned nNewPoints=newPtMap.dataArray.end()->newPointIdx;
	CLOG(logDEBUG) << "newPoints: " << nNewPoints;
	// Allocate space in the new mesh
	newMesh.allocate(nNewPoints);

	/*
	 * Assign all the new points and normals correctly
	 */
	for (unsigned iPoint=0; iPoint < nOgPoints; iPoint++) {

		unsigned newPtIndex=newPtMap.dataArray[iPoint].newPointIdx;
		unsigned ogPtIndex=newPtMap.dataArray[iPoint].pointIdx;
		const float_t * ogCoords = ogMesh.getPointPosition(ogPtIndex);
		const float_t * ogNormal = ogMesh.getNormalVector(ogPtIndex);

		newMesh.setPoint(newPtIndex,ogCoords);
		newMesh.setNormal(newPtIndex,ogNormal);
	}

	/*
	 * Get the triangles
	 */
	unsigned nTriangles = ogMesh.numberOfTriangles();
	std::vector<unsigned> oldToNewMap = newPtMap.oldToNewIdxMap();

	//# pragma omp parallel for
	for (unsigned iTriangle=0; iTriangle < nTriangles;++iTriangle) {
		const unsigned *ogTriangleIdxs=ogMesh.getTriangleIndices(iTriangle);
		unsigned newTriangleIdxs[3];
		newTriangleIdxs[0]=oldToNewMap[ogTriangleIdxs[0]];
		newTriangleIdxs[1]=oldToNewMap[ogTriangleIdxs[1]];
		newTriangleIdxs[2]=oldToNewMap[ogTriangleIdxs[2]];

		newMesh.addTriangle(newTriangleIdxs);
	}

}



#endif /* ALGORITHM_BUILDMESH_H_ */
