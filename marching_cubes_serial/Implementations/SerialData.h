/*
 * SerialData.h
 *
 *  Created on: Aug 17, 2015
 *      Author: sjmunn
 */

#ifndef IMPLEMENTATIONS_SERIALDATA_H_
#define IMPLEMENTATIONS_SERIALDATA_H_

#include"../../common/Algorithm/MarchAlgorithm.h"

template<typename T>
class SerialData : public RuntimeData<T> {
public:
	SerialData();
	virtual ~SerialData();

	void accept(MarchAlgorithm<T> &ma);
};

#endif /* IMPLEMENTATIONS_SERIALDATA_H_ */
