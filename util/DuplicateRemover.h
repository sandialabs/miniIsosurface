/*
 * DuplicateRemover.h
 *
 *  Created on: Jan 25, 2017
 *      Author: dbourge
 */

#ifndef UTIL_DUPLICATEREMOVER_H_
#define UTIL_DUPLICATEREMOVER_H_

#include <array>
#include <vector>
#include <algorithm>

#include "../util/TriangleMesh.h"

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
    std::vector<std::pair<unsigned, unsigned> >& duplicateTracker,
    std::vector<std::array<T, 3> > const&        points,
    std::vector<std::array<T, 3> > const&        normals,
    std::vector<std::array<unsigned, 3> >&       indexTriangles)
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
    unsigned prevEdge = duplicateTracker[0].second;
    duplicateTracker[0].second = 0;
    for(unsigned i = 1; i != duplicateTracker.size(); ++i)
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
    unsigned numPoints = duplicateTracker.back().second + 1;
    std::vector<std::array<T, 3> > mPoints(numPoints);
    std::vector<std::array<T, 3> > mNormals(numPoints);

    std::vector<unsigned> oldToNewMap(duplicateTracker.size());

    for(unsigned i = 0; i != duplicateTracker.size(); ++i)
    {
        unsigned const& oldPointIdx = duplicateTracker[i].first;
        unsigned const& newPointIdx = duplicateTracker[i].second;

        mPoints[newPointIdx] = points[oldPointIdx];
        mNormals[newPointIdx] = normals[oldPointIdx];

        oldToNewMap[oldPointIdx] = newPointIdx;
    }

    // update the triangles using oldToNewMap.
    for(std::array<unsigned, 3>& tri: indexTriangles)
    {
        tri[0] = oldToNewMap[tri[0]];
        tri[1] = oldToNewMap[tri[1]];
        tri[2] = oldToNewMap[tri[2]];
    }

    // mPoints, mNormals and indexTriangles contain all the information
    // needed with respect to this new polygonal mesh, stored in a TriangleMesh.
    return TriangleMesh<T>(mPoints, mNormals, indexTriangles);
}

}

#endif
