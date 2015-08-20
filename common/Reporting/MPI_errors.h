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
class file_not_found: public std::exception {
public:
	file_not_found(const char* const fileName) :
		fileName(fileName) {

		CLOG(logERROR) << "File not found error";
		CLOG(logERROR) << "Given Path: " << fileName;
		CLOG(logERROR) << "In linux, type pwd [Enter] to get your current directory";
		CLOG(logERROR) << "The path should be specified relative to the current directory";
	}

private:
	const char* const fileName;
};



#endif /* REPORTING_MPI_ERRORS_H_ */
