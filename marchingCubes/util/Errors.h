/*
 * Errors.h
 *
 *  Created on: Jul 6, 2015
 *      Author: sjmunn
 *
 * miniIsosurface is distributed under the OSI-approved BSD 3-clause License.
 * See LICENSE.txt for details.
 *
 * Copyright (c) 2017
 * National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under
 * the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
 * certain rights in this software.
 */

#ifndef IO_ERRORS_H_
#define IO_ERRORS_H_

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

class not_enough_memory_allocated: public std::exception {
public:
    not_enough_memory_allocated(const char* const inMessage) :
        message(inMessage) {
            std::cout << "Not enough memory allocated: " << inMessage << std::endl;
    }
private:
    const char* const message;
};

} // util namespace

#endif
