/*
 * Triplet.h
 *
 *  Created on: Jul 14, 2015
 *      Author: sjmunn
 */

/*
 * Triplet Class
 * Purpose: Store triplets of coordinates or indices
 *
 * Usage:   3-dimensional points:
 * 			Triplet<float> point;
 *
 * 			Triangle point index:
 * 			Triplet<unsigned> triangle;
 */

#ifndef DATA_OBJ_TRIPLET_H_
#define DATA_OBJ_TRIPLET_H_

// External includes
#include"../includes.h"

// Reporting Headers
#include"../Reporting/Log.h"

template<typename T>
class Triplet {
public:
	// Constructor/Destructors
	Triplet();
	Triplet(const T, const T, const T);				// Recommended
	Triplet(const T *);
	virtual ~Triplet();

	// Assignment
	void setCoordinates(const T, const T, const T); // Recommended
	void setCoordinates(const T*);					// Risky b/c no bounds check
	void operator=(const T*);
	void operator=(const Triplet&);
	Triplet& operator+=(const T);

	// Comparison
	bool operator==(const Triplet &) const;
	bool operator==(const T*) const;
	bool operator!=(const Triplet &) const;
	bool operator!=(const T*) const;
	bool operator <(const Triplet& rhs) const
	{
		return coordinates[0] < rhs.coordinates[0];
	}

	// Coordinate retrieval
	const T* getCoordinates(void) const;
	const T& operator[](int idx) const { return coordinates[idx]; };
	template<typename Anything> friend std::ostream& operator<<(std::ostream&, const Triplet<Anything>&);

private:
	// Need an array because std::array is a C++11 feature
	T coordinates[3];
};

#endif /* DATA_OBJ_TRIPLET_H_ */
