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
        std::cout << "Given Path: " << fileName << std::endl;
    }

private:
    const char* const fileName;
};

class bad_format: public std::exception {
public:
    bad_format(const char* const inMessage) :
            message(inMessage) {
        std::cout << "bad_format error with message: " << inMessage << std::endl;
    }

private:
    const char* const message;
};

class no_type: public std::exception {
public:
    no_type(const char* const inMessage) :
            message(inMessage) {
        std::cout << "Type not defined error with message: " << inMessage << std::endl;
    }

private:
    const char* const message;
};

} // util namespace

#endif
