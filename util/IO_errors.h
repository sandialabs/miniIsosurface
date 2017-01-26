/*
 * IO_errors.h
 *
 *  Created on: Jul 6, 2015
 *      Author: sjmunn
 */

#ifndef IO_ERRORS_H_
#define IO_ERRORS_H_

//#include"../includes.h"

//#include"Log.h"

#include <iostream>
#include <stdexcept>

using std::size_t;

namespace util {

class file_not_found: public std::exception {
public:
    file_not_found(const char* const inFileName) :
        fileName(inFileName) {

        std::cout << "file_not_found " << inFileName << std::endl;
//      CLOG(logERROR) << "File not found error";
//      CLOG(logERROR) << "Given Path: " << fileName;
//      CLOG(logERROR) << "In linux, type pwd [Enter] to get your current directory";
//      CLOG(logERROR) << "The path should be specified relative to the current directory";
    }

private:
    const char* const fileName;
};

class object_not_sorted: public std::exception {
public:
    object_not_sorted(const char* const inObjectName) :
        objectName(inObjectName) {

        std::cout << "object_not_sorted " << inObjectName << std::endl;
//      CLOG(logERROR) << "The MapReverse object is not sorted";
//      CLOG(logERROR) << "Do not call MapReverse::getNewIndices if the object is not sorted";
    }

private:
    const char* const objectName;
};

class bad_format: public std::exception {
public:
    bad_format(const char* const inMessage) :
            message(inMessage) {
        std::cout << "bad_format " << inMessage << std::endl;
//      CLOG(logERROR) << "Bad format error with message: ";
//      CLOG(logERROR) << message;
    }

private:
    const char* const message;
};

class no_type: public std::exception {
public:
    no_type(const char* const inMessage) :
            message(inMessage) {
        std::cout << "no_type " << inMessage << std::endl;
//      CLOG(logERROR) << "Type not defined error with message: ";
//      CLOG(logERROR) << message;
    }

private:
    const char* const message;
};

class file_too_large: public std::exception {
public:
    file_too_large(const size_t inMessage) :
        message(inMessage) {

        std::cout << "file_too_large " << inMessage << std::endl;
//      CLOG(logERROR) << "Step by Step debugging is active";
//      CLOG(logERROR) << "Only image files up to 1000 points are accepted";
//      CLOG(logERROR) << "The file specified has too many data points";
//      CLOG(logERROR) << "Number of points in file: " << message;
    }

private:
    const size_t message;
};

} // util namespace
#endif
