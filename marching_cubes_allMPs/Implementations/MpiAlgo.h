/*
 * SerialAlgo.h
 *
 *  Created on: Aug 17, 2015
 *      Author: sjmunn
 */

#ifndef IMPLEMENTATIONS_MERGEMPALGO_H_
#define IMPLEMENTATIONS_MERGEMPALGO_H_

// External includes
#include"../../common/includes.h"

// Reporting Headers
#include"../../common/Reporting/YAML_Element.hpp"
#include"../../common/Reporting/YAML_Doc.hpp"
#include"../../common/Reporting/Log.h"
#include"../../common/Reporting/Timer.h"

// Strategy base class
#include"../../common/Algorithm/MarchAlgorithm.h"
#include"../../common/GeneralContext/GeneralContext.h"

#include"../../common/Data_Obj/TriangleMesh.h"

// Algorithm objects
#include"../../common/Algorithm/Ranges.h"

#include"../Algorithm/DuplicateRemover.h"
#include"../Algorithm/buildMesh.h"

// IO Object
#include"../../common/IO/LoadImage3DMPI.h"

//MPI
#include"mpi.h"

template<typename T>
class MpiAlgo : public MarchAlgorithm<T>  {
public:
	MpiAlgo(LoadImage3DMPI<T> &, int, int, Timer *);
	virtual ~MpiAlgo();
	static bool testZeroExtent(unsigned *);

	unsigned numBlocks(const Range oneDRange);
	void march(GeneralContext<T> &);
private:
	// These are intermediate data structures before the final
	// results are sent to GeneralContext
	DuplicateRemover duplicateRemover;
	unsigned grainDim;

	// MPI specific
	LoadImage3DMPI<T> fileHeader;
	int pID,processes;
	Timer *processTimer;
};

#endif /* IMPLEMENTATIONS_MERGEMPALGO_H_ */
