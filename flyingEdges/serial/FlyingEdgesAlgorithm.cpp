/*
 * flyingEdgesAlgorithm.cpp
 *
 *  Created on: Feb 17, 2017
 *      Author: dbourge
 */
#include "FlyingEdgesAlgorithm.h"

#include "../util/MarchingCubesTables.h"

///////////////////////////////////////////////////////////////////////////////
// Pass 1 of the algorithm
///////////////////////////////////////////////////////////////////////////////
void FlyingEdgesAlgorithm::processGridEdges()
{
    for(size_t k = 0; k != nz; ++k) {
    for(size_t j = 0; j != ny; ++j)
    {
        auto curEdgeCases = edgeCases.begin() + (nx-1) * (k*ny + j);
        auto curPointValues = image.getRowIter(j, k);

        gridEdge& curGridEdge = gridEdges[k*ny + j];

        // Calculate xl, xr
        // Fill out edgeCases

        std::array<bool, 2> isGE;
        isGE[0] = (curPointValues[0] >= isoval);
        for(int i = 1; i != nx; ++i)
        {
            isGE[i%2] = (curPointValues[i] >= isoval);

            curEdgeCases[i] = calcCaseEdge(isGE[(i+1)%2], isGE[i%2]);

            // If the edge is cut
            if(curEdgeCases[i] == 1 || curEdgeCases[i] == 2)
            {
                if(curGridEdge.xl == 0)
                    curGridEdge.xl = i;

                curGridEdge.xr = i+1;
            }
        }
    }}
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Pass 2 of the algorithm
///////////////////////////////////////////////////////////////////////////////
void FlyingEdgesAlgorithm::processGridCells()
{
    for(size_t k = 0; k != nz-1; ++k) {
    for(size_t j = 0; j != ny-1; ++j)
    {
        // find adjusted trim values
        size_t xl, xr;
        calcTrimValues(xl, xr, j, k); // xl, xr set in this function

        // comment TODO
        gridEdge& ge0 = gridEdges[k*ny + j];
        gridEdge& ge1 = gridEdges[k*ny + j + 1];
        gridEdge& ge2 = gridEdges[(k+1)*ny + j];
        gridEdge& ge3 = gridEdges[(k+1)*ny + j + 1];

        // comment TODO
        auto const& ec0 = edgeCases.begin() + (nx-1)*(k*ny + j);
        auto const& ec1 = edgeCases.begin() + (nx-1)*(k*ny + j + 1);
        auto const& ec2 = edgeCases.begin() + (nx-1)*((k+1)*ny + j);
        auto const& ec3 = edgeCases.begin() + (nx-1)*((k+1)*ny + j + 1);

        size_t& curTriCounter = *(triCounter.begin() + k*(ny-1) + j);

        auto curCubeCaseIds = cubeCases.begin() + (nx-1)*(k*(ny-1) + j);

        bool isYEnd = (j == ny-2);
        bool isZEnd = (k == nz-2);

        for(size_t i = xl; i != xr; ++i)
        {
            bool isXEnd = (i == nx-2);

            // using edgeCases, compute cubeCases for this cube
            uchar caseId = calcCubeCase(ec0[i], ec1[i], ec2[i], ec3[i]);

            curCubeCaseIds[i] = caseId;

            curTriCounter += util::numTris[caseId];

            const bool* isCut = util::isCut[caseId]; // size 12

            ge0.xstart += isCut[0];
            ge0.ystart += isCut[3];
            ge0.zstart += isCut[8];

            // Note: Each 'gridCell' contains four gridEdges running along it,
            //       ge0, ge1, ge2 and ge3. Each gridCell can access it's own
            //       ge0 but ge1, ge2 and ge3 are owned by other gridCells.
            //       Accessing ge1, ge2 and ge3 leads to a race condition
            //       unless gridCell is along the boundry of the image.

            if(isXEnd)
            {
                ge0.ystart += isCut[1];
                ge0.zstart += isCut[9];
            }
            if(isYEnd)
            {
                ge1.xstart += isCut[2];
                ge1.zstart += isCut[10];
            }
            if(isZEnd)
            {
                ge2.xstart += isCut[4];
                ge2.ystart += isCut[7];
            }

            if(isXEnd and isYEnd)
            {
                ge1.zstart += isCut[11];
            }
            if(isXEnd and isZEnd)
            {
                ge2.ystart += isCut[5];
            }
            if(isYEnd and isZEnd)
            {
                ge3.xstart += isCut[6];
            }
        }
    }}
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Pass 3 of the algorithm
///////////////////////////////////////////////////////////////////////////////
void FlyingEdgesAlgorithm::configureOutputAndAllocate()
{
    // Accumulate triangles into triCounter
    size_t tmp;
    size_t triAccum = 0;
    for(size_t k = 0; k != nz-1; ++k) {
    for(size_t j = 0; j != ny-1; ++j)
    {
        size_t& curTriCounter = triCounter[k*(ny-1)+j];

        tmp = curTriCounter;
        curTriCounter = triAccum;
        triAccum += tmp;
    }}

    // accumulate points
    size_t pointAccum = 0;
    for(size_t k = 0; k != nz; ++k) {
    for(size_t j = 0; j != ny; ++j)
    {
        gridEdge& curGridEdge = gridEdges[k*ny + j];

        tmp = curGridEdge.xstart;
        curGridEdge.xstart = pointAccum;
        pointAccum += tmp;

        tmp = curGridEdge.ystart;
        curGridEdge.ystart = pointAccum;
        pointAccum += tmp;

        tmp = curGridEdge.zstart;
        curGridEdge.zstart = pointAccum;
        pointAccum += tmp;
    }}

    points = std::vector<std::array<float_t, 3> >(pointAccum);
    normals = std::vector<std::array<float_t, 3> >(pointAccum);
    tris = std::vector<std::array<size_t, 3> >(triAccum);
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Pass 4 of the algorithm
///////////////////////////////////////////////////////////////////////////////
void FlyingEdgesAlgorithm::generateOutput()
{
    for(size_t k = 0; k != nz-1; ++k) {
    for(size_t j = 0; j != ny-1; ++j)
    {
        // find adjusted trim values
        size_t xl, xr;
        calcTrimValues(xl, xr, j, k); // xl, xr set in this function

        size_t triIdx = triCounter[k*(ny-1) + j];
        auto curCubeCaseIds = cubeCases.begin() + (nx-1)*(k*(ny-1) + j);

        gridEdge const& ge0 = gridEdges[k*ny + j];
        gridEdge const& ge1 = gridEdges[k*ny + j + 1];
        gridEdge const& ge2 = gridEdges[(k+1)*ny + j];
        gridEdge const& ge3 = gridEdges[(k+1)*ny + j + 1];

        size_t x0counter = 0;
        size_t y0counter = 0;
        size_t z0counter = 0;

        size_t x1counter = 0;
        size_t z1counter = 0;

        size_t x2counter = 0;
        size_t y2counter = 0;

        size_t x3counter = 0;

        bool isYEnd = (j == ny-2);
        bool isZEnd = (k == nz-2);

        for(size_t i = xl; i != xr; ++i)
        {
            bool isXEnd = (i == nx-2);

            uchar caseId = curCubeCaseIds[i]; // is cubeCaseIds even needed?

            if(caseId == 0 || caseId == 255)
            {
                continue;
            }

            const bool* isCut = util::isCut[caseId]; // has 12 elements

            // TODO ineffecient to calculate all of the cube ??
            //      lazy evaluate?
            cube_t        pointCube = image.getPosCube(i, j, k);
            scalarCube_t  isovalCube = image.getValsCube(i, j, k);
            cube_t        gradCube = image.getGradCube(i, j, k);

            // Add Points and normals.
            // Calculate global indices for triangles...

            std::vector<size_t> globalIdxs(12);
            if(isCut[0])
            {
                size_t idx = ge0.xstart + x0counter;
                points[idx] = interpolateOnCube(pointCube, isovalCube, 0);
                normals[idx] = interpolateOnCube(gradCube, isovalCube, 0);
                globalIdxs[0] = idx;
                ++x0counter;
            }

            if(isCut[3])
            {
                size_t idx = ge0.ystart + y0counter;
                points[idx] = interpolateOnCube(pointCube, isovalCube, 3);
                normals[idx] = interpolateOnCube(gradCube, isovalCube, 3);
                globalIdxs[3] = idx;
                ++y0counter;
            }

            if(isCut[8])
            {
                size_t idx = ge0.zstart + z0counter;
                points[idx] = interpolateOnCube(pointCube, isovalCube, 8);
                normals[idx] = interpolateOnCube(gradCube, isovalCube, 8);
                globalIdxs[8] = idx;
                ++z0counter;
            }


            if(isCut[1])
            {
                size_t idx = ge0.ystart + y0counter;
                if(isXEnd)
                {
                    points[idx] = interpolateOnCube(pointCube, isovalCube, 1);
                    normals[idx] = interpolateOnCube(gradCube, isovalCube, 1);
                }
                globalIdxs[1] = idx;
                ++y0counter;
            }

            if(isCut[9])
            {
                size_t idx = ge0.zstart + z0counter;
                if(isXEnd)
                {
                    points[idx] = interpolateOnCube(pointCube, isovalCube, 9);
                    normals[idx] = interpolateOnCube(gradCube, isovalCube, 9);
                }
                globalIdxs[9] = idx;
                ++z0counter;
            }

            if(isCut[2])
            {
                size_t idx = ge1.xstart + x1counter;
                if(isYEnd)
                {
                    points[idx] = interpolateOnCube(pointCube, isovalCube, 2);
                    normals[idx] = interpolateOnCube(gradCube, isovalCube, 2);
                }
                globalIdxs[2] = idx;
                ++x1counter;
            }

            if(isCut[10])
            {
                size_t idx = ge1.zstart + z1counter;
                if(isYEnd)
                {
                    points[idx] = interpolateOnCube(pointCube, isovalCube, 10);
                    normals[idx] = interpolateOnCube(gradCube, isovalCube, 10);
                }
                globalIdxs[10] = idx;
                ++z1counter;
            }

            if(isCut[4])
            {
                size_t idx = ge2.xstart + x2counter;
                if(isZEnd)
                {
                    points[idx] = interpolateOnCube(pointCube, isovalCube, 4);
                    normals[idx] = interpolateOnCube(gradCube, isovalCube, 4);
                }
                globalIdxs[4] = idx;
                ++x2counter;
            }

            if(isCut[7])
            {
                size_t idx = ge2.ystart + y2counter;
                if(isZEnd)
                {
                    points[idx] = interpolateOnCube(pointCube, isovalCube, 7);
                    normals[idx] = interpolateOnCube(gradCube, isovalCube, 7);
                }
                globalIdxs[7] = idx;
                ++y2counter;
            }

            if(isCut[11])
            {
                size_t idx = ge1.zstart + z1counter;
                if(isXEnd and isYEnd)
                {
                    points[idx] = interpolateOnCube(pointCube, isovalCube, 11);
                    normals[idx] = interpolateOnCube(gradCube, isovalCube, 11);
                }
                globalIdxs[11] = idx;
                ++z1counter;
            }

            if(isCut[5])
            {
                size_t idx = ge2.ystart + y2counter;
                if(isXEnd and isZEnd)
                {
                    points[idx] = interpolateOnCube(pointCube, isovalCube, 5);
                    normals[idx] = interpolateOnCube(gradCube, isovalCube, 5);
                }
                globalIdxs[5] = idx;
                ++y2counter;
            }

            if(isCut[6])
            {
                size_t idx = ge3.xstart + x3counter;
                if(isYEnd and isZEnd)
                {
                    points[idx] = interpolateOnCube(pointCube, isovalCube, 6);
                    normals[idx] = interpolateOnCube(gradCube, isovalCube, 6);
                }
                globalIdxs[6] = idx;
                ++x3counter;
            }

            // Add triangles
            const char* caseTri = util::caseTriangles[caseId]; // size 16
            for(int idx = 0; caseTri[idx] != -1; idx += 3)
            {
                tris[triIdx][0] = globalIdxs[caseTri[idx]];
                tris[triIdx][1] = globalIdxs[caseTri[idx+1]];
                tris[triIdx][2] = globalIdxs[caseTri[idx+2]];
                ++triIdx;
            }
        }
    }}
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Don't copy points, normals and tris but move the output into a TrianlgeMesh.
///////////////////////////////////////////////////////////////////////////////
util::TriangleMesh FlyingEdgesAlgorithm::moveOutput()
{
    return util::TriangleMesh(std::move(points),
                              std::move(normals),
                              std::move(tris));
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Private helper functions
///////////////////////////////////////////////////////////////////////////////

inline uchar
FlyingEdgesAlgorithm::calcCaseEdge(
    bool const& prevEdge,
    bool const& currEdge) const
{
    // o -- is greater than or equal to
    // case 0: (i-1) o-----o (i) | (_,j,k)
    // case 1: (i-1) x-----o (i) | (_,j+1,k)
    // case 2: (i-1) o-----x (i) | (_,j,k+1)
    // case 3: (i-1) x-----x (i) | (_,j+1,k+1)
    if(prevEdge && currEdge)
        return 0;
    if(!prevEdge && currEdge)
        return 1;
    if(prevEdge && !currEdge)
        return 2;
    else // !prevEdge && !currEdge
        return 3;
}

inline uchar
FlyingEdgesAlgorithm::calcCubeCase(
    uchar const& ec0, uchar const& ec1,
    uchar const& ec2, uchar const& ec3) const
{
    // ec0 | (_,j,k)
    // ec1 | (_,j+1,k)
    // ec2 | (_,j,k+1)
    // ec3 | (_,j+1,k+1)

    uchar caseId = 0;
    if((ec0 == 0) || (ec0 == 2)) // 0 | (i,j,k)
        caseId |= 1;
    if((ec0 == 0) || (ec0 == 1)) // 1 | (i+1,j,k)
        caseId |= 2;
    if((ec1 == 0) || (ec1 == 1)) // 2 | (i+1,j+1,k)
        caseId |= 4;
    if((ec1 == 0) || (ec1 == 2)) // 3 | (i,j+1,k)
        caseId |= 8;
    if((ec2 == 0) || (ec2 == 2)) // 4 | (i,j,k+1)
        caseId |= 16;
    if((ec2 == 0) || (ec2 == 1)) // 5 | (i+1,j,k+1)
        caseId |= 32;
    if((ec3 == 0) || (ec3 == 1)) // 6 | (i+1,j+1,k+1)
        caseId |= 64;
    if((ec3 == 0) || (ec3 == 2)) // 7 | (i,j+1,k+1)
        caseId |= 128;
    return caseId;
}

inline void
FlyingEdgesAlgorithm::calcTrimValues(
    size_t& xl, size_t& xr,
    size_t const& j, size_t const& k) const
{
    gridEdge const& ge0 = gridEdges[k*ny + j];
    gridEdge const& ge1 = gridEdges[k*ny + j + 1];
    gridEdge const& ge2 = gridEdges[(k+1)*ny + j];
    gridEdge const& ge3 = gridEdges[(k+1)*ny + j + 1];

    using std::min;
    using std::max;
    xl = min(ge0.xl, min(ge1.xl, min(ge2.xl, ge3.xl)));
    xr = max(ge0.xr, max(ge1.xr, max(ge2.xr, ge3.xr)));
}

inline std::array<float_t, 3>
FlyingEdgesAlgorithm::interpolateOnCube(
    cube_t const& pts,
    scalarCube_t const& isovals,
    uchar const& edge) const
{
    uchar i0 = util::edgeVertices[edge][0];
    uchar i1 = util::edgeVertices[edge][1];

    float_t weight = (isoval - isovals[i0]) / (isovals[i1] - isovals[i0]);
    return interpolate(pts[i0], pts[i1], weight);
}

inline std::array<float_t, 3>
FlyingEdgesAlgorithm::interpolate(
    std::array<float_t, 3> const& a,
    std::array<float_t, 3> const& b,
    float_t const& weight) const
{
    std::array<float_t, 3> ret;
    ret[0] = a[0] + (weight * (b[0] - a[0]));
    ret[1] = a[1] + (weight * (b[1] - a[1]));
    ret[2] = a[2] + (weight * (b[2] - a[2]));
    return ret;
}

///////////////////////////////////////////////////////////////////////////////

