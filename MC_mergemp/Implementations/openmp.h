/*
 * openmp.h
 *
 *  Created on: Jul 10, 2015
 *      Author: sjmunn
 */

#ifndef ALGORITHM_OPENMP_H_
#define ALGORITHM_OPENMP_H_

// Common utility headers -------------
// Standard C/C++ library
#include"../../common/includes.h"

// Data Objects
#include"../../common/types.h"

// Reporting Headers
#include"../../common/Reporting/YAML_Element.hpp"
#include"../../common/Reporting/YAML_Doc.hpp"
#include"../../common/Reporting/Log.h"
#include"../../common/Reporting/Timer.h"

// Local implementation headers -------
// Algorithm objects/functions
#include"../Algorithm/Ranges.h"
#include"../Algorithm/marchBlock.h"
#include"../Algorithm/pointCount.h"

namespace openmp {

// Main Marching Cubes Algorithm
void extractIsosurface(const Image3D_t &, float_t, TriangleMesh_t *&, YAML_Doc&);

unsigned numBlocks(const Range);

}

#endif /* ALGORITHM_OPENMP_H_ */
