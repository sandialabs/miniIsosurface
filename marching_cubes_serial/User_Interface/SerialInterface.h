/*
 * MergeMPterface.h
 *
 *  Created on: Aug 3, 2015
 *      Author: sjmunn
 */

#ifndef USER_INTERFACE_SERIALINTERFACE_H_
#define USER_INTERFACE_SERIALINTERFACE_H_

// Common utility headers -------------
// Data Objects
#include"../../common/types.h"

// User Interface
#include"../../common/User_Interface/UI.h"

// Local implementation headers -------
// Marching Cubes Implementation
#include"../Implementations/serial.h"

class SerialInterface: public UI {
public:
	SerialInterface(int argc, char* argv[]) : UI(argc, argv) {};
	virtual ~SerialInterface() {};

	void marchImplemtation(const Image3D_t &, TriangleMesh_t *&, YAML_Doc&);
};

#endif /* USER_INTERFACE_SERIALINTERFACE_H_ */
