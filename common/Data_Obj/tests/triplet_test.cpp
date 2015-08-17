/*
 * triplet_test.cpp
 *
 *  Created on: Jul 14, 2015
 *      Author: sjmunn
 */


/*
 * Tests all passed
 * Ran with
 * g++ -otritest -lboost_unit_test_framework triplet_test.cpp
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE tritest
#include <boost/test/unit_test.hpp>
#include"Triplet.h"
#include"Triplet.cpp" // Bad practice to include .cpp
					  // is necessary for compiling tests though

BOOST_AUTO_TEST_CASE(universeInOrder)
{
		// Assignment array to Triplet
		float coords[]={1.2,1.3,1.4};
		float coordsB[]={2.2,2.3,2.4};
		Triplet<float> ptA(coords);
		ptA=coordsB;
		const float * assignment = ptA.getCoordinates();
		BOOST_CHECK(assignment[0] == coordsB[0]);
		BOOST_CHECK(assignment[1] == coordsB[1]);
		BOOST_CHECK(assignment[2] == coordsB[2]);

		// Assignment of Triplet to Triplet
		Triplet<float> ptB(coords);
		ptA=ptB;
		assignment = ptA.getCoordinates();
		BOOST_CHECK(assignment[0] == coords[0]);
		BOOST_CHECK(assignment[1] == coords[1]);
		BOOST_CHECK(assignment[2] == coords[2]);

		// Comparison of two objects
		BOOST_CHECK(ptA == ptB);

		// Coordinate retrieval
		Triplet<int> pt(1,2,3);
		int coordsC[]={1,2,3};
		const int * othercoords = pt.getCoordinates();
	    BOOST_CHECK(othercoords[0] == coordsC[0]);
	    BOOST_CHECK(othercoords[1] == coordsC[1]);
	    BOOST_CHECK(othercoords[2] == coordsC[2]);

}
