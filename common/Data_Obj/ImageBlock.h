/*
 * ImageBlock.h
 *
 *  Created on: Aug 20, 2015
 *      Author: sjmunn
 */

#ifndef DATA_OBJ_IMAGEBLOCK_H_
#define DATA_OBJ_IMAGEBLOCK_H_

#include "Image3D.h"

template<typename T>
class ImageBlock: public Image3D {
public:
	ImageBlock();
	virtual ~ImageBlock();
};

#endif /* DATA_OBJ_IMAGEBLOCK_H_ */
