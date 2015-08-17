/*
 * serial.h
 *
 *  Created on: Aug 17, 2015
 *      Author: sjmunn
 */

/*
 * Visitor design pattern ---
 * Visited:
 * MCdata		  [parent]			general marching cubes data objects
 * SerialData     [child]			data objects specific to the serial implementation
 * OpenMPData	  [child]			data objects specific to the openMP implementation
 * MergeMPData	  [child]			data objects specific to the openMP with mesh merging implementation
 *
 * Visitors:
 * MarchAlgorithm [parent]			general marching cubes methods
 * SerialAlgo	  [child]			marching cubes methods specific to serial implementation
 * OpenMPAlgo	  [child]			marching cubes methods specific to openMP implementation
 * MergeMPAlgo	  [child]			marching cubes methods specific to openMP with mesh merging implementation
 */

#ifndef SERIAL_H_
#define SERIAL_H_

// External includes
#include"../includes.h"

// Reporting Headers
#include"../Reporting/Log.h"

// Local Includes
#include"../RuntimeData/MCdata.h"

template<typename T>
class serial : public MarchAlgorithm<T> {
public:
	serial();
	virtual ~serial();

	virtual void visit(MCdata<T> *);
};

#endif /* SERIAL_H_ */
