/*
 * LineStream.h
 *
 *  Created on: Jul 6, 2015
 *      Author: sjmunn
 */

#ifndef LINESTREAM_H_
#define LINESTREAM_H_

class LineStream {
public:
	LineStream(std::istream &in);
	std::istream& stream();

	void readline();

private:
	std::istream &in;
	std::stringstream sstream;
	std::string line;
};

inline LineStream::LineStream(std::istream &input) :
		in(input) {
}

inline std::istream& LineStream::stream() {
	return this->sstream;
}

inline void LineStream::readline() {
	std::getline(this->in, this->line);
	this->sstream.clear();
	this->sstream.str(this->line);
	this->sstream.seekg(0);
}

#endif /* LINESTREAM_H_ */
