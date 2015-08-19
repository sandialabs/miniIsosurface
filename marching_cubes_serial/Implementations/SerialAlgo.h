/*
 * SerialAlgo.h
 *
 *  Created on: Aug 17, 2015
 *      Author: sjmunn
 */

#ifndef IMPLEMENTATIONS_SERIALALGO_H_
#define IMPLEMENTATIONS_SERIALALGO_H_

// External includes
#include"../../common/includes.h"

// Reporting Headers
#include"../../common/Reporting/YAML_Element.hpp"
#include"../../common/Reporting/YAML_Doc.hpp"
#include"../../common/Reporting/Log.h"

// Strategy base class
#include"../../common/Algorithm/MarchAlgorithm.h"
#include"../../common/GeneralContext/GeneralContext.h"

// Algorithm objects
#include"../../common/Algorithm/Ranges.h"

template<typename T>
class SerialAlgo : public MarchAlgorithm<T>  {
public:
	SerialAlgo();
	virtual ~SerialAlgo();

	void march(GeneralContext<T> &);
};

#endif /* IMPLEMENTATIONS_SERIALALGO_H_ */
