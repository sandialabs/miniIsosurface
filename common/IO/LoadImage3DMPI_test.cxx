/*
 * LoadImage3DMPI_test.cxx
 *
 *  Created on: Aug 20, 2015
 *      Author: sjmunn
 */

/*
 * LoadImage3DMPI_test.cpp
 *
 *  Created on: Jul 16, 2015
 *      Author: sjmunn
 */

/*
 * Tests all passed
 * Ran with
 * export LD_LIBRARY_PATH=/usr/lib64
 * g++ -std=c++0x -oiotest -lboost_unit_test_framework LoadImage3DMPI_test.cxx
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE tritest
#include <boost/test/unit_test.hpp>

#include"LoadImage3DMPI.h"
#include"LoadImage3DMPI.cpp" // Bad practice to include .cpp
// is necessary for compiling tests though

#include"../Data_Obj/Image3D.h"
#include"../Data_Obj/Image3D.cpp"
#include"../Reporting/YAML_Doc.hpp"
#include"../Reporting/YAML_Doc.cpp"
#include"../Reporting/YAML_Element.hpp"
#include"../Reporting/YAML_Element.cpp"

BOOST_AUTO_TEST_CASE(readingHeader) {
	LoadImage3DMPI<float_t> headerInfo;
	headerInfo.loadHeader("../../Data/ctmini.vtk");
	const unsigned *dims = headerInfo.getVolumeDimensions();
	BOOST_CHECK(dims[0] == 3);
	BOOST_CHECK(dims[1] == 3);
	BOOST_CHECK(dims[2] == 3);
	BOOST_CHECK(headerInfo.getnVolumePoints() == 27);

}

BOOST_AUTO_TEST_CASE(readingFirstBlock) {
	BOOST_TEST_CHECKPOINT("Getting header info");
	LoadImage3DMPI<float_t> headerInfo;
	headerInfo.loadHeader("../../Data/ctmini.vtk");
	const unsigned *dims = headerInfo.getVolumeDimensions();

	BOOST_TEST_CHECKPOINT("Initializing new reader object");
	LoadImage3DMPI<float_t> firstBlock(headerInfo);
	BOOST_TEST_CHECKPOINT("firstBlock is initialized to binary data in the file");
	unsigned extent[6] = { 0, 0, 0, 0, 0, 0 };
	firstBlock.setBlockExtent(extent);
	Image3D<float_t> imageVol;
	BOOST_TEST_CHECKPOINT("Read the image block");
	firstBlock.readBlockData(imageVol);

	// get the point from imageVol
	const float_t *buffer;
	buffer = imageVol.getData();
	for (int iPoint=0;iPoint<8;++iPoint) {
		std::cout << "buffer val is: " << buffer[iPoint] << std::endl;
	}
}

BOOST_AUTO_TEST_CASE(readingAnotherBlock) {
	BOOST_TEST_CHECKPOINT("Getting header info");
	LoadImage3DMPI<float_t> headerInfo;
	headerInfo.loadHeader("../../Data/ctmini.vtk");
	const unsigned *dims = headerInfo.getVolumeDimensions();

	BOOST_TEST_CHECKPOINT("Initializing new reader object");
	LoadImage3DMPI<float_t> firstBlock(headerInfo);
	BOOST_TEST_CHECKPOINT("firstBlock is initialized to binary data in the file");
	unsigned extent[6] = { 1, 1, 1, 1, 1, 1 };
	firstBlock.setBlockExtent(extent);
	Image3D<float_t> imageVol;
	BOOST_TEST_CHECKPOINT("Read the image block");
	firstBlock.readBlockData(imageVol);

	// get the point from imageVol
	const float_t *buffer;
	buffer = imageVol.getData();
	for (int iPoint=0;iPoint<8;++iPoint) {
		std::cout << "buffer val is: " << buffer[iPoint] << std::endl;
	}
}
