/*
 * serial_test.cpp
 *
 *  Created on: Aug 17, 2015
 *      Author: sjmunn
 */

/*
 * Testing in progress
 * Ran with
 * export LD_LIBRARY_PATH=/usr/lib64
 * g++ -oserialtest -lboost_unit_test_framework serial_test.cxx
 * -- when all Data_Objs are included
 * g++ -std=c++0x -lrt -oserialtest -lboost_unit_test_framework serial_test.cxx
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE tritest
#include <boost/test/unit_test.hpp>

#include "serial.h"

// Because there's no makefile
#include "serial.cxx"
#include "../RuntimeData/MCdata.cxx"
#include "MarchAlgorithm.cpp"
#include "../RuntimeData/RuntimeData.cpp"
// Data Object definitions
#include"../Data_Obj/Image3D.cpp"
#include"../Data_Obj/TriangleMesh.cpp"
#include"../Data_Obj/Triplet.cpp"
#include"../Data_Obj/MapReverse.cxx"
#include"../Data_Obj/EdgeIndexer.cpp"
// Reporting Headers
#include"../Reporting/YAML_Element.cpp"
#include"../Reporting/YAML_Doc.cpp"
#include"../Reporting/Timer.cpp"

BOOST_AUTO_TEST_CASE(universeInOrder)
{
		MCdata<float_t> data;
		serial<float_t> serialObj;
		data.accept(serialObj);
}
