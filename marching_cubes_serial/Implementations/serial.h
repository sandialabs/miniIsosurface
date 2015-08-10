/*
 * serial.h
 *
 *  Created on: Jul 7, 2015
 *      Author: sjmunn
 */

#ifndef ALGORITHM_SERIAL_H_
#define ALGORITHM_SERIAL_H_

// Common utility headers -------------
// Standard C/C++ library
#include"../../utils/includes.h"

// Data Objects
#include"../../utils/types.h"

// Reporting Headers
#include"../../utils/Reporting/YAML_Element.hpp"
#include"../../utils/Reporting/YAML_Doc.hpp"
#include"../../utils/Reporting/Log.h"
#include"../../utils/Reporting/Timer.h"

// Local implementation headers -------
// Algorithm objects/functions
#include"../Algorithm/Ranges.h"
#include"../Algorithm/marchBlock.h"

namespace serial {

// Main Marching Cubes Algorithm
void extractIsosurface(const Image3D_t &, float_t, TriangleMesh_t *&, YAML_Doc&);

}

#endif /* ALGORITHM_SERIAL_H_ */
