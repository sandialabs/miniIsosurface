/*
 * Ranges.h
 *
 *  Created on: Jul 10, 2015
 *      Author: sjmunn
 */

#ifndef ALGORITHM_RANGES_H_
#define ALGORITHM_RANGES_H_

class Range {
public:
	Range() :
			from(0), to(0), gsize(1) {
	}
	Range(unsigned inFrom, unsigned inTo, unsigned inGsize = 1) :
			from(inFrom), to(inTo), gsize(inGsize) {
	}

	unsigned begin() const {
		return from;
	}

	unsigned end() const {
		return to;
	}

	unsigned grain() const {
		return gsize;
	}

private:
	unsigned from, to, gsize;
};

class Range2D {
public:
	Range2D() :
			r(0, 0), c(0, 0) {
	}
	Range2D(unsigned rfrom, unsigned rto, unsigned rgrain, unsigned cfrom,
			unsigned cto, unsigned cgrain) :
			r(rfrom, rto, rgrain), c(cfrom, cto, cgrain) {
	}
	Range2D(unsigned rfrom, unsigned rto, unsigned cfrom, unsigned cto) :
			r(rfrom, rto), c(cfrom, cto) {
	}

	const Range& rows() const {
		return r;
	}

	const Range& cols() const {
		return c;
	}

private:
	Range r, c;
};

class Range3D {
public:
	Range3D() :
			p(0, 0), r(0, 0), c(0, 0) {
	}
	Range3D(unsigned pfrom, unsigned pto, unsigned pgrain, unsigned rfrom,
			unsigned rto, unsigned rgrain, unsigned cfrom, unsigned cto,
			unsigned cgrain) :
			p(pfrom, pto, pgrain), r(rfrom, rto, rgrain), c(cfrom, cto, cgrain) {
	}
	Range3D(unsigned pfrom, unsigned pto, unsigned rfrom, unsigned rto,
			unsigned cfrom, unsigned cto) :
			p(pfrom, pto), r(rfrom, rto), c(cfrom, cto) {
	}

	const Range& pages() const {
		return p;
	}

	const Range& rows() const {
		return r;
	}

	const Range& cols() const {
		return c;
	}

	void extent(unsigned * inExtent) const {
		inExtent[0] = this->c.begin();
		inExtent[1] = this->c.end() - 1;
		inExtent[2] = this->r.begin();
		inExtent[3] = this->r.end() - 1;
		inExtent[4] = this->p.begin();
		inExtent[5] = this->p.end() - 1;
	}

private:
	Range p, r, c;
};

#endif /* ALGORITHM_RANGES_H_ */
