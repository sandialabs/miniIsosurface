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

template<typename T>
class Image3D  {
public:
	Image3D();
	virtual ~Image3D();

	void setDimension(unsigned, unsigned, unsigned);
	void setSpacing(T, T, T);
	void setOrigin(T, T, T);

	virtual void allocate();
	T *getData();

	const unsigned* getDimension() const;
	unsigned getNumberOfPoints() const;
	const T* getSpacing() const;
	const T* getOrigin() const;
	const T* getData() const;

	void report(YAML_Doc &) const;

private:
	// Prevent object copying
	Image3D(const Image3D&); // no implementation
	Image3D& operator=(const Image3D&); // no implementation

	unsigned dim[3];
	unsigned npoints;
	T spacing[3];
	T origin[3];
	T *data;
};

#endif /* DATA_OBJ_IMAGE3D_H_ */
