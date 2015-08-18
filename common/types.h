/*
 * types.h
 *
 *  Created on: Jul 10, 2015
 *      Author: sjmunn
 */

#ifndef TYPES_H_
#define TYPES_H_

// Standard library
#include <map>

// Data Object definitions
#include"./Data_Obj/Image3D.h"
#include"./Data_Obj/TriangleMesh.h"
#include"./Data_Obj/Triplet.h"
#include"./Data_Obj/MapReverse.h"
#include"./Algorithm/MarchAlgorithm.h"
#include"./Algorithm/EdgeIndexer.h"


// Data Object names
typedef Image3D<float_t> Image3D_t;
typedef TriangleMesh<float_t> TriangleMesh_t;
typedef EdgeIndexer<float_t> EdgeIndexer_t;

typedef Triplet<float_t> PositionVector_t;
typedef Triplet<unsigned> IndexTriplet_t;

typedef std::unordered_map<unsigned,unsigned> PointMap_t;

//typedef MarchAlgorithm<float_t> MarchAlgorithm_t;
//typedef GeneralContext<float_t> GeneralContext_t;

#endif /* TYPES_H_ */
