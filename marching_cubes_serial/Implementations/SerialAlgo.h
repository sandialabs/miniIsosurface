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
#include"../../common/Reporting/Log.h"

// Local Includes
#include"SerialData.h"

template<typename T>
class SerialAlgo : public MarchAlgorithm<T>  {
public:
	SerialAlgo();
	virtual ~SerialAlgo();

	virtual void visit(MCdata<T> *);
};

#endif /* IMPLEMENTATIONS_SERIALALGO_H_ */
