/*
 * MPterface.h
 *
 *  Created on: Aug 3, 2015
 *      Author: sjmunn
 */

#ifndef USER_INTERFACE_MPTERFACE_H_
#define USER_INTERFACE_MPTERFACE_H_

// Common utility headers -------------
// Data Objects
#include"../../utils/types.h"

// User Interface
#include"../../utils/User_Interface/UI.h"

// Local implementation headers -------
// Marching Cubes Implementation
#include"../Implementations/openmp.h"

class MPterface: public UI {
public:
	MPterface(int argc, char* argv[]) : UI(argc, argv) {};
	virtual ~MPterface() {};

	void marchImplemtation(const Image3D_t &, TriangleMesh_t *&, YAML_Doc&);
};

#endif /* USER_INTERFACE_MPTERFACE_H_ */
