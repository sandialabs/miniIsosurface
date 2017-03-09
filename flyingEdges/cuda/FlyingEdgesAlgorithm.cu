/*
 * flyingEdgesAlgorithm.cpp
 *
 *  Created on: Feb 17, 2017
 *      Author: dbourge
 */
#include "FlyingEdgesAlgorithm.h"

#include "CudaMarchingCubesTables.h"

#include <numeric>

#include <iostream> // TODO

// TODO make sure pointValues stored in const memory

///////////////////////////////////////////////////////////////////////////////
// Pass 1 of the algorithm
///////////////////////////////////////////////////////////////////////////////

__device__
uchar calcCaseEdge(
    bool const& prevEdge,
    bool const& currEdge)
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

__global__
void pass1gpu_edgeCases(
    scalar_t* pointValues,
    scalar_t isoval,
    int nx, int ny,
    uchar* edgeCases)
{
    // (nx-1, ny, nz) > comes as (nx-1, ny*nz)
    // Each row has several blocks
    // Each thread is one point

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y % ny;
    int k = blockIdx.y / ny;

    __shared__ bool isGE[FE_BLOCK_WIDTH_PLUS_ONE];

    if(i < nx)
        isGE[threadIdx.x] = pointValues[k*nx*ny + j*nx + i] >= isoval;

    if(threadIdx.x == 0 && i < nx-1)
    {
        isGE[blockDim.x] = pointValues[k*nx*ny + j*nx + i + blockDim.x] >= isoval;
    }

    __syncthreads();

    if(i < nx-1)
    {
        uchar caseEdge = calcCaseEdge(isGE[threadIdx.x], isGE[threadIdx.x + 1]);
        edgeCases[k*(nx-1)*ny + j*(nx-1) + i] = caseEdge;
    }
}

__global__
void pass1gpu_trim(
    int nx, int ny, int nz,                    // input
    uchar* edgeCases,                          // input
    FlyingEdgesAlgorithm::gridEdge* gridEdges) // output
{
    // (1, ny, nz) > comes as (ny, nz)

    int j = blockIdx.x * blockDim.x + threadIdx.x;
    int k = blockIdx.y * blockDim.y + threadIdx.y;

    if(j >= ny || k >= nz)
        return;

    size_t xl = nx;
    size_t xr = 0;

    uchar* curEdgeCases = edgeCases + k*(nx-1)*ny + j*(nx-1);

    for(int i = 0; i != nx-1; ++i)
    {
        if(curEdgeCases[i] == 1 || curEdgeCases[i] == 2)
        {
            if(xl == nx)
                xl = i;
            xr = i+1;
        }
    }

    gridEdges[k*ny + j].xl = xl;
    gridEdges[k*ny + j].xr = xr;
}

void FlyingEdgesAlgorithm::pass1()
{
    int tx = FE_BLOCK_WIDTH;
    uint3 gridDim = make_uint3(((nx-1) + tx - 1) / tx, ny*nz, 1);
    uint3 blockDim = make_uint3(tx, 1, 1);
    pass1gpu_edgeCases<<<gridDim, blockDim>>>(
        pointValues,
        isoval,
        nx, ny,
        edgeCases);

    int ty = FE_BLOCK_WIDTH_Y;
    int tz = FE_BLOCK_WIDTH_Z;
    gridDim = make_uint3((ny + ty - 1) / ty, (nz + tz - 1) / tz, 1);
    blockDim = make_uint3(ty, tz, 1);

    pass1gpu_trim<<<gridDim, blockDim>>>(
        nx, ny, nz,
        edgeCases,
        gridEdges);

    cudaDeviceSynchronize();

    /////////////////////////////////////
    // WHAT IS GOING ON?..It works now //  TODO
    /////////////////////////////////////

    int numGE = nz*ny;
    gridEdge* hostGEs = (gridEdge*)malloc(numGE*sizeof(gridEdge));
    cudaMemcpy(hostGEs, gridEdges, numGE*sizeof(gridEdge),
               cudaMemcpyDeviceToHost);

    size_t countL = 0;
    size_t countR = 0;
    for(int idx = 0; idx != numGE; ++idx)
    {
        countL += hostGEs[idx].xl;
        countR += hostGEs[idx].xr;
    }

    std::cout << "xl, xr: " << countL << ", " << countR << std::endl;

    free(hostGEs);
}

/*
void FlyingEdgesAlgorithm::pass1()
{
    // For each (j, k):
    //  - for each edge i along fixed (j, k) gridEdge, fill edgeCases with
    //    cut information.
    //  - find the locations for computational trimming, xl and xr
    for(size_t k = 0; k != nz; ++k) {
    for(size_t j = 0; j != ny; ++j)
    {
        auto curEdgeCases = edgeCases.begin() + (nx-1) * (k*ny + j);
        auto curPointValues = image.getRowIter(j, k);

        gridEdge& curGridEdge = gridEdges[k*ny + j];

        std::array<bool, 2> isGE;
        isGE[0] = (curPointValues[0] >= isoval);
        for(int i = 1; i != nx; ++i)
        {
            isGE[i%2] = (curPointValues[i] >= isoval);

            curEdgeCases[i-1] = calcCaseEdge(isGE[(i+1)%2], isGE[i%2]);

            // If the edge is cut
            if(curEdgeCases[i-1] == 1 || curEdgeCases[i-1] == 2)
            {
                if(curGridEdge.xl == 0)
                    curGridEdge.xl == i-1;

                curGridEdge.xr = i;
            }
        }
    }}
}

///////////////////////////////////////////////////////////////////////////////
*/
///////////////////////////////////////////////////////////////////////////////
// Pass 2 of the algorithm
///////////////////////////////////////////////////////////////////////////////

__device__
void calcTrimValues(
    int& xl, int& xr,
    FlyingEdgesAlgorithm::gridEdge const& ge0,
    FlyingEdgesAlgorithm::gridEdge const& ge1,
    FlyingEdgesAlgorithm::gridEdge const& ge2,
    FlyingEdgesAlgorithm::gridEdge const& ge3)
{
    xl = min(ge0.xl, min(ge1.xl, min(ge2.xl, ge3.xl)));
    xr = max(ge0.xr, max(ge1.xr, max(ge2.xr, ge3.xr)));
}

__device__
uchar calcCubeCase(
    uchar const& ec0, uchar const& ec1,
    uchar const& ec2, uchar const& ec3)
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

__global__
void pass2gpu_cubeCases(
    int nx, int ny, int nz,
    uchar* edgeCases,
    FlyingEdgesAlgorithm::gridEdge* gridEdges,
    int* triCounter,
    uchar* cubeCases)
{
    // (1, ny-1, nz-1) > comes as (ny-1, nz-1)
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    int k = blockIdx.y * blockDim.y + threadIdx.y;

    if(j >= ny-1 || k >= nz-1)
        return;

    FlyingEdgesAlgorithm::gridEdge& ge0 = gridEdges[k*ny + j];
    FlyingEdgesAlgorithm::gridEdge& ge1 = gridEdges[k*ny + j + 1];
    FlyingEdgesAlgorithm::gridEdge& ge2 = gridEdges[(k+1)*ny + j];
    FlyingEdgesAlgorithm::gridEdge& ge3 = gridEdges[(k+1)*ny + j + 1];

    uchar* ec0 = edgeCases + k*ny*(nx-1) + j*(nx-1); // (nx-1)*(k*ny + j);
    uchar* ec1 = edgeCases + k*ny*(nx-1) + (j+1)*(nx-1); // (nx-1)*(k*ny + j + 1);
    uchar* ec2 = edgeCases + (k+1)*ny*(nx-1) + j*(nx-1); // (nx-1)*((k+1)*ny + j);
    uchar* ec3 = edgeCases + (k+1)*ny*(nx-1) + (j+1)*(nx-1); // (nx-1)*((k+1)*ny + j + 1);

    int xl, xr;
    calcTrimValues(xl, xr, ge0, ge1, ge2, ge3);

    int triCount = 0;
    uchar* curCubeCases = cubeCases + k*(nx-1)*(ny-1) + j*(nx-1);

    int xstart = 0;
    int ystart = 0;
    int zstart = 0; // TODO don't set initial values in gridEdge Constructor;

    const bool* isCut;
    for(int i = xl; i != xr; ++i) // What happens here on a gpu?
                                  // I imagine it takes the max xr-xl of all blocks
    {
        // TODO why is this needed?
        if (i >= nx-1)
            return;

        uchar caseId = calcCubeCase(ec0[i], ec1[i], ec2[i], ec3[i]);

        curCubeCases[i] = caseId; // THIS LINE BREAKS EVERYTHING

        // Can't imagine this would do anything on a gpu unless all threads
        // on a block evaluated to the same value.
        if(caseId == 0 || caseId == 255)
        {
            continue;
        }

        triCount += cuda_util::numTris[caseId];
        isCut = cuda_util::isCut[caseId]; // if xr == nx-1, then xr-1 is cut
                                          // so this will be set

        xstart += isCut[0];
        ystart += isCut[3];
        zstart += isCut[8];
    }

    triCounter[k*(ny-1) + j] = triCount;

    if(xr == nx-1)
    {
        // isCut was set at i = xr-1
        ystart += isCut[1];
        zstart += isCut[9];
    }

    ge0.xstart = xstart;
    ge0.ystart = ystart;
    ge0.zstart = zstart;
}

__global__
void pass2gpu_ghost_xz(
    int nx, int ny, int nz,
    uchar* edgeCases,
    FlyingEdgesAlgorithm::gridEdge* gridEdges)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;

    //if(k >= nz-1)
    //    return;
    if(k >= nz) // This function will deal with gridEdge at (_, ny-1, nz-1)
        return;

    bool isCorner = k == nz-1;

    int j = ny-1;

    FlyingEdgesAlgorithm::gridEdge& ge0 = gridEdges[k*ny + j];
    // If isCorner, this is just bogus.
    FlyingEdgesAlgorithm::gridEdge& ge1 = gridEdges[(1-isCorner)*(k+1)*ny + j];

    uchar* ec0 = edgeCases + k*ny*(nx-1) + j*(nx-1);
    // If isCorner, this is just bogus
    uchar* ec1 = edgeCases + (1-isCorner)*(k+1)*ny*(nx-1) + j*(nx-1);

    int xl = min(ge0.xl, nx*isCorner + (1-isCorner)*ge1.xl);
    int xr = max(ge0.xr, (1-isCorner)*ge1.xr);

    if(xl >= xr)
        return;

    int xstart = 0;
    int zstart = 0; // TODO don't set initial values in gridEdge Constructor;

    uchar c0;
    uchar c1;

    for(int i = xl; i != xr; ++i)
    {
        c0 = ec0[i];
        c1 = ec1[i];

        // see if the edges are cut
        xstart += (c0 == 1 || c0 == 2);

        // bogus if isCorner
        zstart += ( (c0 == 0 && c1 == 1) || (c0 == 0 && c1 == 3) ||
                    (c0 == 1 && c1 == 2) || (c0 == 2 && c1 == 3) );
    }

    if(xr == nx-1)
    {
        // bogus if isCorner
        zstart += ( (c0 == 0 && c1 == 2) || (c0 == 0 && c1 == 3) ||
                    (c0 == 1 && c1 == 2) || (c0 == 1 && c1 == 3) );
    }

    ge0.xstart = xstart;
    ge0.ystart = 0;
    ge0.zstart = zstart*(1-isCorner);
}

__global__
void pass2gpu_ghost_xy(
    int nx, int ny, int nz,
    uchar* edgeCases,
    FlyingEdgesAlgorithm::gridEdge* gridEdges)
{
    int j = blockIdx.x * blockDim.x + threadIdx.x;

    if(j >= ny-1)
        return;

    int k = nz-1;

    FlyingEdgesAlgorithm::gridEdge& ge0 = gridEdges[k*ny + j];
    FlyingEdgesAlgorithm::gridEdge& ge1 = gridEdges[k*ny + j + 1];

    uchar* ec0 = edgeCases + k*ny*(nx-1) + j*(nx-1);
    uchar* ec1 = edgeCases + k*ny*(nx-1) + (j+1)*(nx-1);

    int xl = min(ge0.xl, ge1.xl);
    int xr = max(ge0.xr, ge1.xr);

    if(xl >= xr)
        return;

    int xstart = 0;
    int ystart = 0; // TODO don't set initial values in gridEdge Constructor;

    uchar c0;
    uchar c1;

    for(int i = xl; i != xr; ++i)
    {
        c0 = ec0[i];
        c1 = ec1[i];

        // see if the edges are cut
        xstart += (c0 == 1 || c0 == 2);
        ystart += ( (c0 == 0 && c1 == 1) || (c0 == 0 && c1 == 3) ||
                    (c0 == 1 && c1 == 2) || (c0 == 2 && c1 == 3) );
    }

    if(xr == nx-1)
    {
        ystart += ( (c0 == 0 && c1 == 2) || (c0 == 0 && c1 == 3) ||
                    (c0 == 1 && c1 == 2) || (c0 == 1 && c1 == 3) );
    }

    ge0.xstart = xstart;
    ge0.ystart = ystart;
    ge0.zstart = 0;
}

// TOO SLOW! done in xz ghost function
//__global__
//void pass2gpu_ghost_xyz(
//    int nx, int ny, int nz,
//    uchar* edgeCases,
//    FlyingEdgesAlgorithm::gridEdge* gridEdges)
//{
//    int j = ny-1;
//    int k = nz-1;
//
//    FlyingEdgesAlgorithm::gridEdge& ge = gridEdges[k*ny + j];
//    uchar* ec = edgeCases + k*ny*(nx-1) + j*(nx-1);
//
//    int xl = ge.xl;
//    int xr = ge.xr;
//
//    int xstart = 0;
//
//    uchar c;
//
//    for(int i = xl; i != xr; ++i)
//    {
//        c = ec[i];
//        xstart += (c == 1 || c == 2);
//    }
//
//    ge.xstart = xstart;
//    ge.ystart = 0;
//    ge.zstart = 0;
//}

void FlyingEdgesAlgorithm::pass2()
{
    // pass2 calculates
    //   1) cubeCases for each block ray
    //   2) triCount for each block ray
    //   3) edgeRay count

    // 1st kernel: Calculate the 0, 1, 2 edge ray, cube cases, tricount
    // 2nd kernel: Calculate lost edges

    int ty = FE_BLOCK_WIDTH_Y;
    int tz = FE_BLOCK_WIDTH_Z;
    uint3 gridDim = make_uint3(((ny-1) + ty - 1) / ty, ((nz-1) + tz - 1) / tz, 1);
    uint3 blockDim = make_uint3(ty, tz, 1);

    pass2gpu_cubeCases<<<gridDim, blockDim>>>(
        nx, ny, nz,
        edgeCases,
        gridEdges,   // modified
        triCounter,  // modified
        cubeCases);  // modified
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    size_t sz = (nx-1)*(ny-1)*(nz-1)*sizeof(uchar);
    cudaDeviceSynchronize();
    uchar* hostCubeCases = (uchar*)malloc(sz);
    cudaMemcpy(hostCubeCases, cubeCases,
               sz, cudaMemcpyDeviceToHost);

    int count = 0;
    // TODO hostCubeCases is not the same every time.
    for(int i = 0; i != (nx-1)*(ny-1)*(nz-1); ++i)
    {
        if(hostCubeCases[i] != 0 && hostCubeCases[i] != 255)
            count += hostCubeCases[i];
    }
    std::cout << "Count cube cases " << count << std::endl;

    free(hostCubeCases);
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    // TODO these can be launched and executed independently of each other
    int bw = FE_BLOCK_WIDTH;

    // Making sure that the xz face takes care of the (_, ny-1, nz-1) gridEdge
    // BE CAREFUL. xz takes care of corner. don't use (nz-1)
    pass2gpu_ghost_xz<<<(nz + bw - 1) / bw, bw>>>(
        nx, ny, nz,
        edgeCases,
        gridEdges);
    pass2gpu_ghost_xy<<<((ny-1) + bw - 1) / bw, bw>>>(
        nx, ny, nz,
        edgeCases,
        gridEdges);

    cudaDeviceSynchronize();

// This is prohibitively slow so pass2gpu_ghost_xz covers it now
//    pass2gpu_ghost_xyz<<<1, 1, 3>>>(
//        nx, ny, nz,
//        edgeCases,
//        gridEdges);
}

/*
void FlyingEdgesAlgorithm::pass2()
{
    // For each (j, k):
    //  - for each cube (i, j, k) calculate caseId and number of gridEdge cuts
    //    in the x, y and z direction.
    for(size_t k = 0; k != nz-1; ++k) {
    for(size_t j = 0; j != ny-1; ++j)
    {
        // find adjusted trim values
        size_t xl, xr;
        calcTrimValues(xl, xr, j, k); // xl, xr set in this function

        // ge0 is owned by this (i, j, k). ge1, ge2 and ge3 are only used for
        // boundary cells.
        gridEdge& ge0 = gridEdges[k*ny + j];
        gridEdge& ge1 = gridEdges[k*ny + j + 1];
        gridEdge& ge2 = gridEdges[(k+1)*ny + j];
        gridEdge& ge3 = gridEdges[(k+1)*ny + j + 1];

        // ec0, ec1, ec2 and ec3 were set in pass 2. They are used
        // to calculate the cell caseId.
        auto const& ec0 = edgeCases.begin() + (nx-1)*(k*ny + j);
        auto const& ec1 = edgeCases.begin() + (nx-1)*(k*ny + j + 1);
        auto const& ec2 = edgeCases.begin() + (nx-1)*((k+1)*ny + j);
        auto const& ec3 = edgeCases.begin() + (nx-1)*((k+1)*ny + j + 1);

        // Count the number of triangles along this row of cubes.
        size_t& curTriCounter = *(triCounter.begin() + k*(ny-1) + j);

        auto curCubeCaseIds = cubeCases.begin() + (nx-1)*(k*(ny-1) + j);

        bool isYEnd = (j == ny-2);
        bool isZEnd = (k == nz-2);

        for(size_t i = xl; i != xr; ++i)
        {
            bool isXEnd = (i == nx-2);

            // using edgeCases from pass 2, compute cubeCases for this cube
            uchar caseId = calcCubeCase(ec0[i], ec1[i], ec2[i], ec3[i]);

            curCubeCaseIds[i] = caseId;

            // If the cube has no triangles through it
            if(caseId == 0 || caseId == 255)
            {
                continue;
            }

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
            //
            //       To really make sense of the indices, it helps to draw
            //       out the following picture of a cube with the appropriate
            //       labels:
            //         v0 is at (i,   j,   k)
            //         v1       (i+1, j,   k)
            //         v2       (i+1, j+1, k)
            //         v3       (i,   j+1, k)
            //         v4       (i,   j,   k+1)
            //         v5       (i+1, j,   k+1)
            //         v6       (i+1, j+1, k+1)
            //         v7       (i,   j+1, k+1)
            //         e0  connects v0 to v1 and is parallel to the x-axis
            //         e1           v1    v2                        y
            //         e2           v2    v3                        x
            //         e3           v0    v3                        y
            //         e4           v4    v5                        x
            //         e5           v5    v6                        y
            //         e6           v6    v7                        x
            //         e7           v4    v7                        y
            //         e8           v0    v4                        z
            //         e9           v1    v5                        z
            //         e10          v3    v7                        z
            //         e11          v2    v6                        z

            // Handle cubes along the edge of the image
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
*/
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Pass 3 of the algorithm
///////////////////////////////////////////////////////////////////////////////

__global__
void pass3gpu_blockAccum(
    int nx, int ny, int nz, // which are needed TODO?
    int* triCounter,
    FlyingEdgesAlgorithm::gridEdge* gridEdges,
    int* blockAccum)
{
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    int k = blockIdx.y * blockDim.y + threadIdx.y;

    if (j == 0 && k == 0);
    {
        // TODO get rid of this exp
        blockAccum[0] = 191230;
        blockAccum[1] = 192340;
        blockAccum[2] = 193450;
        blockAccum[3] = 194560;
    }
    return;

    // step 1: accumulate individual y thread
    // step 2: calc block sum
    // step 3: __syncthreads
    // step 4: add to individual y thread


    if(k >= nz)
        return;

    __shared__ int accum[4*FE_BLOCK_WIDTH];

    int tmp;
    int accumX   = 0;
    int accumY   = 0;
    int accumZ   = 0;
    int accumTri = 0;
    for(int j = 0; j != ny; ++j)
    {
        FlyingEdgesAlgorithm::gridEdge& ge = gridEdges[k*ny + j];

        tmp = ge.xstart;
        ge.xstart = accumX;
        accumX += tmp;

        tmp = ge.ystart;
        ge.ystart = accumY;
        accumY += tmp;

        tmp = ge.zstart;
        ge.zstart = accumZ;
        accumZ += tmp;
    }

    if(k < nz-1)
    {
        for(int j = 0; j != ny-1; ++j)
        {
            int& curTriCount = triCounter[k*(ny-1) + j];

            tmp = curTriCount;
            curTriCount = accumTri;
            accumTri += tmp;
        }
    }

    accum[4*threadIdx.z + 0] = accumX;
    accum[4*threadIdx.z + 1] = accumY;
    accum[4*threadIdx.z + 2] = accumZ;
    accum[4*threadIdx.z + 3] = accumTri;

    __syncthreads();

    if(threadIdx.z == 0) // agh!
    {
        for(int idx = 1; idx != blockDim.z; ++idx)
        {
            accum[4*idx + 0] += accum[4*(idx-1) + 0];
            accum[4*idx + 1] += accum[4*(idx-1) + 1];
            accum[4*idx + 2] += accum[4*(idx-1) + 2];
            accum[4*idx + 3] += accum[4*(idx-1) + 3];
        }

        // answer for global accumulation
        blockAccum[4*blockIdx.z + 0] = accum[4*(blockDim.z-1) + 0];
        blockAccum[4*blockIdx.z + 1] = accum[4*(blockDim.z-1) + 1];
        blockAccum[4*blockIdx.z + 2] = accum[4*(blockDim.z-1) + 2];
        blockAccum[4*blockIdx.z + 3] = accum[4*(blockDim.z-1) + 3];
    }

    __syncthreads();

    if(threadIdx.z == 0)
        return;

    bool isEndK = k == nz-1;
    for(int j = 1; j != ny-1; ++j)
    {
        FlyingEdgesAlgorithm::gridEdge& ge = gridEdges[k*ny + j];

        ge.xstart += accum[4*(threadIdx.z-1) + 0];
        ge.ystart += accum[4*(threadIdx.z-1) + 1];
        ge.zstart += accum[4*(threadIdx.z-1) + 2];

        // put z stuff here..
        if(!isEndK)
            triCounter[k*(ny-1) + j] = accum[4*(threadIdx.z-1) + 3];
    }

    FlyingEdgesAlgorithm::gridEdge& ge = gridEdges[k*ny + (ny-1)];
    ge.xstart += accum[4*(threadIdx.z-1) + 0];
    ge.ystart += accum[4*(threadIdx.z-1) + 1];
    ge.zstart += accum[4*(threadIdx.z-1) + 2];
}

__global__ // TODO can split up along j here easy enough.
void pass3gpu_gridAccum(
    int nx, int ny, int nz, // which are needed TODO?
    int* triCounter,
    FlyingEdgesAlgorithm::gridEdge* gridEdges,
    int* blockAccum) // used as input here
{
    // not adding to the first block!
    //
    // add to individual y threads
    int k = (blockIdx.z + 1)*blockDim.z + threadIdx.z;

    if (k >= nz)
        return;

    int addX   = blockAccum[4*blockIdx.z + 0];
    int addY   = blockAccum[4*blockIdx.z + 1];
    int addZ   = blockAccum[4*blockIdx.z + 2];
    int addTri = blockAccum[4*blockIdx.z + 3];

    for(int j = 0; j != ny; ++j)
    {
        FlyingEdgesAlgorithm::gridEdge& ge = gridEdges[k*ny + j];
        ge.xstart += addX;
        ge.ystart += addY;
        ge.zstart += addZ;
    }

    if(k >= nz-1)
        return;

    for(int j = 0; j != ny-1; ++j)
    {
        triCounter[k*(ny-1) + j] += addTri;
    }
}

// Can make prettier?
void FlyingEdgesAlgorithm::pass3()
{
    // Split the z axis
    // Kernel 1: calculate the accum values on block sync
    //           then accum individual values
    // Use that info accum each block (except the first one)
    // Kernel 2: just add values to individual threads
    int tz = FE_BLOCK_WIDTH;

    int numBlocks = (nz + tz - 1) / tz;

    // there are four because: xstart, ystart, zstart, triaccum
    int sizeBlocks = 4 * numBlocks * sizeof(int);

    uint3 gridDim = make_uint3(1, numBlocks, 1); // TODO FIGURE OUT HOW THREADS PER BLOCK
                                                 //      STUFF WORKS with at
                                                 //      3rd dimension..
                                                 //
                                                 //      Blocks can have 3
                                                 //      dimensions
                                                 //
                                                 //      Grids can only have 2
                                                 //      dimensions!
    uint3 blockDim = make_uint3(1, tz, 1);

    int* hostBlockAccum = (int*)malloc(sizeBlocks);

    int* deviceBlockAccum;
    cudaMalloc(&deviceBlockAccum, sizeBlocks);

    std::cout << gridDim.x << ", " << gridDim.y << std::endl;
    std::cout << blockDim.x << ", " << blockDim.y << std::endl;

    // Accumulate values locally

    pass3gpu_blockAccum<<<gridDim, blockDim>>>(
        nx, ny, nz,
        triCounter,
        gridEdges,
        deviceBlockAccum);

    cudaDeviceSynchronize();

    cudaMemcpy(hostBlockAccum, deviceBlockAccum,
               sizeBlocks, cudaMemcpyDeviceToHost);


    //////////////////////////////////////////////////////////////
    // WHAT IS GOING ON HERE?
    //////////////////////////////////////////////////////////////

    std::cout  << "SANITY CHECK " << hostBlockAccum[0] << ", " << hostBlockAccum[1] << ", " << hostBlockAccum[2]  << ", " << hostBlockAccum[3] << std::endl;
//    if(err != cudaSuccess)
//        std::cout << "AGHHHHHHHHHHHHHHHH" << std::endl;
//    else
//        std::cout << __LINE__ << std::endl;
//    for(int idx = 0; idx != numBlocks; ++idx)
//    {
//        std::cout << hostBlockAccum[4*idx + 3] << std::endl;
//    }

    if(numBlocks != 1)
    {

        // std::partial_sum(2 2 3 4  3  2  2 ) TODO not using it get rid of header
        // goes to         (2 4 7 11 14 16 18)
        // std::partial_sum(hostBlockAccum, hostBlockAccum + numBlocks, hostBlockAccum);

        for(int i = 4; i != 4*numBlocks; i += 4)
        {
            hostBlockAccum[i+0] += hostBlockAccum[i-4];
            hostBlockAccum[i+1] += hostBlockAccum[i-3];
            hostBlockAccum[i+2] += hostBlockAccum[i-2];
            hostBlockAccum[i+3] += hostBlockAccum[i-1];
        }
        // note: the last values in hostBlockAccum should contain total counts

        // The first block is done so it is ignored
        // and the last info in BlockAccum isn't needed (its the total counts)
        cudaMemcpy(deviceBlockAccum, hostBlockAccum,
                   sizeBlocks - 4 * sizeof(int), cudaMemcpyHostToDevice);

        // Accumulate values from other blocks
        gridDim = make_uint3(1, 1, numBlocks - 1);
        pass3gpu_gridAccum<<<gridDim, blockDim>>>(
            nx, ny, nz,
            triCounter,
            gridEdges,
            deviceBlockAccum);
    }

    // Allocate memory for points, normals and tris
    outputAllocated = true;
    size_t numPoints = hostBlockAccum[4*(numBlocks-1) + 0] +
                       hostBlockAccum[4*(numBlocks-1) + 1] +
                       hostBlockAccum[4*(numBlocks-1) + 2];
    size_t numTris   = hostBlockAccum[4*(numBlocks-1) + 3];

    cudaMalloc(&points,  3*sizeof(scalar_t)*numPoints);
    cudaMalloc(&normals, 3*sizeof(scalar_t)*numPoints);
    cudaMalloc(&tris, 3*sizeof(int)*numTris);

    std::cout << "PASS3 " << numPoints << " " << numTris << std::endl;

    // free memory used in this function
    free(hostBlockAccum);
    cudaFree(deviceBlockAccum);

    cudaDeviceSynchronize();

}
/*
void FlyingEdgesAlgorithm::pass3()
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

    // accumulate points, filling out starting locations of each gridEdge
    // in the process.
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

    points = std::vector<std::array<scalar_t, 3> >(pointAccum);
    normals = std::vector<std::array<scalar_t, 3> >(pointAccum);
    tris = std::vector<std::array<size_t, 3> >(triAccum);
}
*/
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Pass 4 of the algorithm
///////////////////////////////////////////////////////////////////////////////
void FlyingEdgesAlgorithm::pass4()
{
/* Copy of pass 2, should be similar, just different kernels
 *
    int ty = FE_BLOCK_WIDTH_Y;
    int tz = FE_BLOCK_WIDTH_Z;
    uint3 gridDim = make_uint3(1, ((ny-1) + ty - 1) / ty, ((nz-1) + tz - 1) / tz);
    uint3 blockDim = make_uint3(1, ty, tz);

    pass2gpu_cubeCases<<<gridDim, blockDim>>>(
        nx, ny, nz,
        edgeCases,
        gridEdges,   // modified
        triCounter,  // modified
        cubeCases);  // modified

    // TODO these can be launched and executed independently of each other
    int bw = FE_BLOCK_WIDTH;

    // Making sure that the xz face takes care of the (_, ny-1, nz-1) gridEdge
    // BE CAREFUL. xz takes care of corner. don't use (nz-1)
    pass2gpu_ghost_xz<<<(nz + bw - 1) / bw, bw>>>(
        nx, ny, nz,
        edgeCases,
        gridEdges);
    pass2gpu_ghost_xy<<<((ny-1) + bw - 1) / bw, bw>>>(
        nx, ny, nz,
        edgeCases,
        gridEdges);
*/
}
/*
void FlyingEdgesAlgorithm::pass4()
{
    // For each (j, k):
    //  - For each cube at i, fill out points, normals and triangles owned by
    //    the cube. Each cube is in charge of filling out e0, e3 and e8. Only
    //    in edge cases does it also fill out other edges.
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

            uchar caseId = curCubeCaseIds[i];

            if(caseId == 0 || caseId == 255)
            {
                continue;
            }

            const bool* isCut = util::isCut[caseId]; // has 12 elements

            // Most of the information contained in pointCube, isovalCube
            // and gradCube will be used--but not necessarily all. It has
            // not been tested whether or not obtaining only the information
            // needed will provide a significant speedup--but
            // most likely not.
            cube_t        pointCube = image.getPosCube(i, j, k);
            scalarCube_t  isovalCube = image.getValsCube(i, j, k);
            cube_t        gradCube = image.getGradCube(i, j, k);

            // Add Points and normals.
            // Calculate global indices for triangles
            std::array<size_t, 12> globalIdxs;
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

            // Note:
            //   e1, e5, e9 and e11 will be visited in the next iteration
            //   when they are e3, e7, e8 and 10 respectively. So don't
            //   increment their counters. When the cube is an edge cube,
            //   their counters don't need to be incremented because they
            //   won't be used agin.

            // Manage boundary cases if needed. Otherwise just update
            // globalIdx.
            if(isCut[1])
            {
                size_t idx = ge0.ystart + y0counter;
                if(isXEnd)
                {
                    points[idx] = interpolateOnCube(pointCube, isovalCube, 1);
                    normals[idx] = interpolateOnCube(gradCube, isovalCube, 1);
                    // y0counter counter doesn't need to be incremented
                    // because it won't be used again.
                }
                globalIdxs[1] = idx;
            }

            if(isCut[9])
            {
                size_t idx = ge0.zstart + z0counter;
                if(isXEnd)
                {
                    points[idx] = interpolateOnCube(pointCube, isovalCube, 9);
                    normals[idx] = interpolateOnCube(gradCube, isovalCube, 9);
                    // z0counter doesn't need to in incremented.
                }
                globalIdxs[9] = idx;
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
                    // z1counter does not need to be incremented.
                }
                globalIdxs[11] = idx;
            }

            if(isCut[5])
            {
                size_t idx = ge2.ystart + y2counter;
                if(isXEnd and isZEnd)
                {
                    points[idx] = interpolateOnCube(pointCube, isovalCube, 5);
                    normals[idx] = interpolateOnCube(gradCube, isovalCube, 5);
                    // y2 counter does not need to be incremented.
                }
                globalIdxs[5] = idx;
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

inline std::array<scalar_t, 3>
FlyingEdgesAlgorithm::interpolateOnCube(
    cube_t const& pts,
    scalarCube_t const& isovals,
    uchar const& edge) const
{
    uchar i0 = util::edgeVertices[edge][0];
    uchar i1 = util::edgeVertices[edge][1];

    scalar_t weight = (isoval - isovals[i0]) / (isovals[i1] - isovals[i0]);
    return interpolate(pts[i0], pts[i1], weight);
}

inline std::array<scalar_t, 3>
FlyingEdgesAlgorithm::interpolate(
    std::array<scalar_t, 3> const& a,
    std::array<scalar_t, 3> const& b,
    scalar_t const& weight) const
{
    std::array<scalar_t, 3> ret;
    ret[0] = a[0] + (weight * (b[0] - a[0]));
    ret[1] = a[1] + (weight * (b[1] - a[1]));
    ret[2] = a[2] + (weight * (b[2] - a[2]));
    return ret;
}

///////////////////////////////////////////////////////////////////////////////
*/
