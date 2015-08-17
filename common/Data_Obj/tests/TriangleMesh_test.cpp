/*
 * TriangleMesh_test.cpp
 *
 *  Created on: Jul 16, 2015
 *      Author: sjmunn
 */


/*
 * Tests all passed
 * Ran with
 * g++ -omeshtest -lboost_unit_test_framework TriangleMesh_test.cpp
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE tritest
#include <boost/test/unit_test.hpp>
#include"Triplet.h"
#include"Triplet.cpp" // Bad practice to include .cpp
					  // is necessary for compiling tests though
#include"TriangleMesh.h"
#include"TriangleMesh.cpp"

BOOST_AUTO_TEST_CASE(universeInOrder)
{
		TriangleMesh<float_t> TheMesh;
		float ptA[]={1.1,1.2,1.3};
		float ptB[]={2.1,2.2,2.3};
		float norm[]={3.1,3.2,3.3};

		TheMesh.addPoint(ptA);
		TheMesh.addPoint(ptB[0],ptB[1],ptB[2]);
		TheMesh.addNormal(norm);

		unsigned idx[]={0,1,2};
		TheMesh.addTriangle(idx);

		Triplet<float_t> pointA(ptA), pointB(ptB);
		Triplet<float_t> normA(norm);
		Triplet<unsigned> triangle(idx);
		std::cout << "The mesh indices," << std::endl;
		std::cout << TheMesh.getTriangleIndices(0)[0] << std::endl;
		std::cout << TheMesh.getTriangleIndices(0)[1] << std::endl;
		std::cout << TheMesh.getTriangleIndices(0)[2] << std::endl;
		BOOST_CHECK(pointA == TheMesh.getPointPosition(0));
		BOOST_CHECK(pointB == TheMesh.getPointPosition(1));
		BOOST_CHECK(normA == TheMesh.getNormalVector(0));
		BOOST_CHECK(triangle == TheMesh.getTriangleIndices(0));


		BOOST_CHECK(TheMesh);
}
