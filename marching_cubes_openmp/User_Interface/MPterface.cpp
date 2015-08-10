/*
 * MergeMPterface.cpp
 *
 *  Created on: Aug 3, 2015
 *      Author: sjmunn
 */

#include "MPterface.h"

void MPterface::marchImplemtation(const Image3D_t & vol, TriangleMesh_t *&mesh,
		YAML_Doc &doc) {
	CLOG(logDEBUG) << "Called implementation";
	CLOG(logDEBUG) << UI::getIsoVal();
	openmp::extractIsosurface(vol, UI::getIsoVal(), mesh, doc);
}
