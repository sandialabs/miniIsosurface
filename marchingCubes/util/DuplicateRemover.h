/*
 * DuplicateRemover.h
 *
 *  Created on: Jan 25, 2017
 *      Author: dbourge
 *
 * miniIsosurface is distributed under the OSI-approved BSD 3-clause License.
 * See LICENSE.txt for details.
 *
 * Copyright (c) 2017
 * National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under
 * the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
 * certain rights in this software.
 */

#ifndef UTIL_DUPLICATEREMOVER_H_
#define UTIL_DUPLICATEREMOVER_H_

#include <array>
#include <vector>
#include <algorithm>

#include "../util/TriangleMesh.h"

using std::size_t;

namespace util {

struct PairCompareSecond
{
    template <typename A, typename B>
    bool operator()(std::pair<A, B> const& lhs, std::pair<A, B> const& rhs)
    {
        return lhs.second < rhs.second;
    }
};

template <typename T>
TriangleMesh<T>
duplicateRemover(
    std::vector<std::pair<size_t, size_t> >& duplicateTracker,
    std::vector<std::array<T, 3> > const&        points,
    std::vector<std::array<T, 3> > const&        normals,
    std::vector<std::array<size_t, 3> >&       indexTriangles)
{
    // Sort by global edge indices
    std::sort(
        duplicateTracker.begin(), duplicateTracker.end(),
        PairCompareSecond());

    // If two subsequent global edge indices are the same, then that
    // point is a duplicate. This block of code rewrites over the
    // global edge indices with new point indices.
    // So if duplicate Tracker has (0, 1), (3, 2), (2, 2), (4, 20) at
    // first, this block of code will leave duplicate tracker with
    // (0, 0), (3, 1), (2, 1), (4, 2)
    size_t prevEdge = duplicateTracker[0].second;
    duplicateTracker[0].second = 0;
    for(size_t i = 1; i != duplicateTracker.size(); ++i)
    {
        if(prevEdge == duplicateTracker[i].second)
        {
            duplicateTracker[i].second = duplicateTracker[i-1].second;
        }
        else
        {
            prevEdge = duplicateTracker[i].second;
            duplicateTracker[i].second = duplicateTracker[i-1].second + 1;
        }
    }

    // Rewrite to new points and normals vector without any
    // duplicates.
    // And create a map from old indices to new indices.
    size_t numPoints = duplicateTracker.back().second + 1;
    std::vector<std::array<T, 3> > mPoints(numPoints);
    std::vector<std::array<T, 3> > mNormals(numPoints);

    std::vector<size_t> oldToNewMap(duplicateTracker.size());

    for(size_t i = 0; i != duplicateTracker.size(); ++i)
    {
        size_t const& oldPointIdx = duplicateTracker[i].first;
        size_t const& newPointIdx = duplicateTracker[i].second;

        mPoints[newPointIdx] = points[oldPointIdx];
        mNormals[newPointIdx] = normals[oldPointIdx];

        oldToNewMap[oldPointIdx] = newPointIdx;
    }

    // update the triangles using oldToNewMap.
    for(std::array<size_t, 3>& tri: indexTriangles)
    {
        tri[0] = oldToNewMap[tri[0]];
        tri[1] = oldToNewMap[tri[1]];
        tri[2] = oldToNewMap[tri[2]];
    }

    // mPoints, mNormals and indexTriangles contain all the information
    // needed with respect to this new polygonal mesh, stored in a TriangleMesh.
    return TriangleMesh<T>(mPoints, mNormals, indexTriangles);
}

} // util namespace

#endif
