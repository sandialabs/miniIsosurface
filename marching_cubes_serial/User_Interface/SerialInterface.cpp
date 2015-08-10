/*
 * MergeMPterface.cpp
 *
 *  Created on: Aug 3, 2015
 *      Author: sjmunn
 */

#include "SerialInterface.h"

void SerialInterface::marchImplemtation(const Image3D_t & vol, TriangleMesh_t *&mesh,
		YAML_Doc &doc) {
	serial::extractIsosurface(vol, UI::getIsoVal(), mesh, doc);
}
