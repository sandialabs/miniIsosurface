/*
 * FlyingEdges_Config.h
 *
 * miniIsosurface is distributed under the OSI-approved BSD 3-clause License.
 * See LICENSE.txt for details.
 *
 * Copyright (c) 2017
 * National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under
 * the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
 * certain rights in this software.
 */

#ifndef FLYINGEDGES_CONFIG_H_
#define FLYINGEDGES_CONFIG_H_

#include <array>

using std::size_t;

using scalar_t = float;

using uchar = unsigned char;

using cube_t = std::array<std::array<scalar_t, 3>, 8>;
using scalarCube_t = std::array<scalar_t, 8>;

#define FE_BLOCK_WIDTH 512
#define FE_BLOCK_WIDTH_PLUS_ONE 513

// FE_BLOCK_WIDTH = FE_BLOCK_WIDTH_Y * FE_BLOCK_WIDTH_Z
#define FE_BLOCK_WIDTH_Y 16
#define FE_BLOCK_WIDTH_Z 32

#endif
