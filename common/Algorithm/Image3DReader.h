/*
 * Image3DReader.h
 *
 *  Created on: Aug 21, 2015
 *      Author: sjmunn
 */

#ifndef IMAGE3DREADER_H_
#define IMAGE3DREADER_H_

#include"../includes.h"

// Reporting Headers
#include"../Reporting/YAML_Element.hpp"
#include"../Reporting/YAML_Doc.hpp"
#include"../Reporting/Log.h"

// Data Objects
#include"../Data_Obj/Image3D.h"

template<typename T>
class Image3DReader {
public:
	Image3DReader(Image3D<T>*);
	virtual ~Image3DReader();

	// Read the Image3D object
	void setImage3DOutputBuffers(unsigned, unsigned,unsigned);
	void getVertexValues(T *,unsigned,unsigned);

	void getValsForGradient(T (& x)[3][2], T (& run)[3], const unsigned, const unsigned, const unsigned) const;

	Image3D<T> *imageData;
private:
	// Prevent object copying (would cost too much memory
	Image3DReader(const Image3DReader&); // no implementation
	Image3DReader& operator=(const Image3DReader&); // no implementation

	/* ===== Reading image object ===============
	 * Collection of buffers to allow the Image3D
	 * object to read itself efficiently.
	 ==========================================*/
	unsigned bufferIdx;
	/*
	 * 4 buffers improve cache efficiency
	 * this speeds up run time significantly
	 */
	const T *X1buffer;
	const T *X2buffer;
	const T *X3buffer;
	const T *X4buffer;

	bool isMPIimage;
};

#endif /* IMAGE3DREADER_H_ */
