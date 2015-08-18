/*
 * GeneralContext.h
 *
 *  Created on: Aug 17, 2015
 *      Author: sjmunn
 */

#ifndef GENERALCONTEXT_H_
#define GENERALCONTEXT_H_

#include"../includes.h"

// Data Object definitions
#include"../Data_Obj/Image3D.h"
#include"../Data_Obj/TriangleMesh.h"
#include"../Data_Obj/Triplet.h"
#include"../Data_Obj/MapReverse.h"

#include"../Algorithm/MarchAlgorithm.h"

template<typename T>
class GeneralContext {
public:
	GeneralContext();
	GeneralContext(MarchAlgorithm<T> *);
	virtual ~GeneralContext();

	void setAlgorithm(MarchAlgorithm<T> *);

	// Data writing methods
	Image3D<T> * imgAdr(void) { return &imageIn;};

	void march(void);
public:
	// General/minimal marching cubes runtime data
	Image3D<T> imageIn;
	unsigned ext[6];
	T isoval;
	TriangleMesh<T> mesh;

	// Initialize console log and YAML
	YAML_Doc doc;
private:
	// Prevent object copying
	GeneralContext(const GeneralContext&); // no implementation
	GeneralContext& operator=(const GeneralContext&); // no implementation

	MarchAlgorithm<T> *_strategy;
};

#endif /* GENERALCONTEXT_H_ */
