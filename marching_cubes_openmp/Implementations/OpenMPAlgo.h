/*
 * SerialAlgo.h
 *
 *  Created on: Aug 17, 2015
 *      Author: sjmunn
 */

#ifndef IMPLEMENTATIONS_OPENMPALGO_H_
#define IMPLEMENTATIONS_OPENMPALGO_H_

// External includes
#include"../../common/includes.h"

// Reporting Headers
#include"../../common/Reporting/YAML_Element.hpp"
#include"../../common/Reporting/YAML_Doc.hpp"
#include"../../common/Reporting/Log.h"

// Strategy base class
#include"../../common/Algorithm/MarchAlgorithm.h"
#include"../../common/GeneralContext/GeneralContext.h"

// Local Includes
#include"../Algorithm/Ranges.h"

template<typename T>
class OpenMPAlgo : public MarchAlgorithm<T>  {
public:
	OpenMPAlgo();
	virtual ~OpenMPAlgo();

	unsigned numBlocks(const Range oneDRange);
	void march(GeneralContext<T> &);
};

#endif /* IMPLEMENTATIONS_OPENMPALGO_H_ */
