/*
 * Image3D.h
 *
 *  Created on: Jul 6, 2015
 *      Author: sjmunn
 */

#ifndef DATA_OBJ_IMAGE3D_H_
#define DATA_OBJ_IMAGE3D_H_

#include"../includes.h"

// Reporting Headers
#include"../Reporting/YAML_Element.hpp"
#include"../Reporting/YAML_Doc.hpp"
#include"../Reporting/Log.h"

template<typename U> class Image3DReader;

template<typename T>
class Image3D  {
public:
	Image3D();
	virtual ~Image3D();

	void setDimension(unsigned, unsigned, unsigned);
	void setMPIorigin(unsigned, unsigned, unsigned);
	void setSpacing(T, T, T);
	void setOrigin(T, T, T);
	void setToMPIdataBlock(void);

	void allocate();
	T *getData();

	const unsigned* getDimension() const;
	unsigned getNumberOfPoints() const;
	const T* getSpacing() const;
	const T* getOrigin() const;
	const T* getData() const;

//	void report(YAML_Doc &) const;

	template<typename anything> friend class Image3DReader;
//	void setImage3DOutputBuffers(const unsigned, const unsigned,const unsigned);
//	void getVertexValues(T *,unsigned,unsigned);
//
//	void getValsForGradient(T (& x)[3][2], const unsigned, const unsigned, const unsigned) const;

private:
	// Prevent object copying (would cost too much memory
	Image3D(const Image3D&); // no implementation
	Image3D& operator=(const Image3D&); // no implementation

	// ===== Image data =========================
	unsigned dim[3];
	unsigned sliceSize;
	unsigned npoints;
	T spacing[3];
	T origin[3];
	T *data;

	bool isMPIdataBlock;
	unsigned MPIorigin[3];

//	/* ===== Reading image object ===============
//	 * Collection of buffers to allow the Image3D
//	 * object to read itself efficiently.
//	 ==========================================*/
//	unsigned bufferIdx;
//	/*
//	 * 4 buffers improve cache efficiency
//	 * this speeds up run time significantly
//	 */
//	const T *X1buffer;
//	const T *X2buffer;
//	const T *X3buffer;
//	const T *X4buffer;
//
//	/* ===== Reading image object ===============
//	 * Collection of functions to allow the Image3D
//	 * object to read itself efficiently.
//	 ==========================================*/
};

#endif /* DATA_OBJ_IMAGE3D_H_ */
