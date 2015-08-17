/*
 * UI.h
 *
 *  Created on: Jul 27, 2015
 *      Author: sjmunn
 */

#ifndef USER_INTERFACE_UI_H_
#define USER_INTERFACE_UI_H_
// Standard C/C++ library
#include"../includes.h"

// Reporting
#include"../Reporting/RunTime_errors.h"
#include"../Reporting/YAML_Element.hpp"
#include"../Reporting/YAML_Doc.hpp"
#include"../Reporting/Log.h"

template<typename T>
class UI {
public:
	UI(int argc, char* argv[]);
	virtual ~UI() {} ;

	//void marchImplemtation(const Image3D_t &, TriangleMesh_t *&, YAML_Doc&);
	const char * getFile(void) const;
	const char * outFile(void) const;
	T getIsoVal(void) const;

private:
	static void checkArgs(int);
	void checkDebugLevel(void) const;

	T isoval;
	char * filePath;
	char * outFilePath;
};

#endif /* USER_INTERFACE_UI_H_ */
