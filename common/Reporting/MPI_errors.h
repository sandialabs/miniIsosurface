/*
 * MPI_erros.h
 *
 *  Created on: Aug 20, 2015
 *      Author: sjmunn
 */

#ifndef REPORTING_MPI_ERRORS_H_
#define REPORTING_MPI_ERRORS_H_

#include"../includes.h"

#include"Log.h"

/*
 * LoadImage3DMPI errors
 */
class zero_dimensions: public std::exception {
public:
	zero_dimensions(const char* const message) :
		message(message) {

		CLOG(logERROR) << "When initializing LoadImage3DMPI object A from another LoadImage3DMPI object B,";
		CLOG(logERROR) << "argument B was found to have zero dimensions.";
		CLOG(logERROR) << "Either B did not yet read its file's header";
		CLOG(logERROR) << "or the file has not actual image volume data.";
	}

private:
	const char* const message;
};

class block_extent_not_set: public std::exception {
public:
	block_extent_not_set(const char* const message) :
		message(message) {

		CLOG(logERROR) << "Cannot call readBlockData() before setting the block extent via setBlockExtent()";
	}

private:
	const char* const message;
};

class impossible_extent: public std::exception {
public:
	impossible_extent(const char* const message) :
		message(message) {

		CLOG(logERROR) << "In LoadImage3D setBlockExtent() member function,";
		CLOG(logERROR) << message;
	}

private:
	const char* const message;
};



#endif /* REPORTING_MPI_ERRORS_H_ */
