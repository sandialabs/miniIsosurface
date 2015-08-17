/*
 * SerialAlgo.h
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

#ifndef IMPLEMENTATIONS_SERIALALGO_H_
#define IMPLEMENTATIONS_SERIALALGO_H_

// External includes
#include"../../common/includes.h"

// Reporting Headers
#include"../../common/Reporting/YAML_Element.hpp"
#include"../../common/Reporting/YAML_Doc.hpp"
#include"../../common/Reporting/Log.h"
#include"../../common/Reporting/Timer.h"


// Local Includes
#include"./SerialData.h"
#include"../Algorithm/Ranges.h"

template<typename T>
class SerialAlgo : public MarchAlgorithm<T>  {
public:
	SerialAlgo();
	virtual ~SerialAlgo();

	void visit(SerialData<T> &);
};

#endif /* IMPLEMENTATIONS_SERIALALGO_H_ */
