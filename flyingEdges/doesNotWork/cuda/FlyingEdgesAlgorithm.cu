/*
 * flyingEdgesAlgorithm.cpp
 *
 *  Created on: Feb 17, 2017
 *      Author: dbourge
 */
#include "FlyingEdgesAlgorithm.h"

#include "CudaMarchingCubesTables.h"

#include <numeric>

#include <algorithm> // TODO
#include <iostream> // TODO

#define MAX_X_GRID 65535 // says larger but doesnt work if larger..
#define MAX_Y_GRID 65535
#define MAX_Z_GRID 65535
#define MAX_X_BLOCK 1024
#define MAX_Y_BLOCK 1024
#define MAX_Z_BLOCK 64
#define MAX_THREAD_PER_BLOCK 1024

#define DEBUG true

// TODO figure out how to handle errors

bool validKernelSize(uint3 const& gridDim, uint3 const& blockDim)
{
    if(gridDim.x > MAX_X_GRID)
        return false;
    if(gridDim.y > MAX_Y_GRID)
        return false;
    if(gridDim.z > MAX_Z_GRID)
        return false;

    if(blockDim.x > MAX_X_BLOCK)
        return false;
    if(blockDim.y > MAX_Y_BLOCK)
        return false;
    if(blockDim.z > MAX_Z_BLOCK)
        return false;

    if(blockDim.x * blockDim.y * blockDim.z > MAX_THREAD_PER_BLOCK)
        return false;

    return true;
}

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
    // Each row has several blocks
    // Each thread is one point

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y;
    int k = blockIdx.z;

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
    uint3 gridDim = make_uint3(((nx-1) + tx - 1) / tx, ny, nz);
    uint3 blockDim = make_uint3(tx, 1, 1);

    if(!validKernelSize(gridDim, blockDim))
        std::cout << "GAHHHHHHHHHHHH GP2 ENGINE " << __LINE__ << std::endl; // TODO

    pass1gpu_edgeCases<<<gridDim, blockDim>>>(
        pointValues,
        isoval,
        nx, ny,
        edgeCases);

    int ty = FE_BLOCK_WIDTH_Y;
    int tz = FE_BLOCK_WIDTH_Z;
    gridDim = make_uint3((ny + ty - 1) / ty, (nz + tz - 1) / tz, 1);
    blockDim = make_uint3(ty, tz, 1);

    if(!validKernelSize(gridDim, blockDim))
        std::cout << "GAHHHHHHHHHHHH GP2 ENGINE " << __LINE__ << std::endl; // TODO

    pass1gpu_trim<<<gridDim, blockDim>>>(
        nx, ny, nz,
        edgeCases,
        gridEdges);

    cudaDeviceSynchronize();

    if(DEBUG)
    {
        int numGE = nz*ny;
        gridEdge* hostGEs = (gridEdge*)malloc(numGE*sizeof(gridEdge));
        cudaMemcpy(hostGEs, gridEdges, numGE*sizeof(gridEdge),
                   cudaMemcpyDeviceToHost);

        int numCubes=(nx-1)*ny*nz;
        size_t count = 0;
        uchar* hoseEdgeCases = (uchar*)malloc(numCubes*sizeof(uchar));
        cudaMemcpy(hoseEdgeCases, edgeCases, numCubes*sizeof(uchar),
                   cudaMemcpyDeviceToHost);
        for(int idx = 0; idx != numCubes; ++idx)
        {
            uchar const& val = hoseEdgeCases[idx];
            count += val;
        }
        std::cout << "Edgecase counter: " << count << std::endl;
        free(hoseEdgeCases);


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
}

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

    if(xl > xr)
        xl = xr;
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
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    int k = blockIdx.y * blockDim.y + threadIdx.y;

    if(j >= ny-1 || k >= nz-1)
        return;

    FlyingEdgesAlgorithm::gridEdge& ge0 = gridEdges[k*ny + j];
    FlyingEdgesAlgorithm::gridEdge& ge1 = gridEdges[k*ny + j + 1];
    FlyingEdgesAlgorithm::gridEdge& ge2 = gridEdges[(k+1)*ny + j];
    FlyingEdgesAlgorithm::gridEdge& ge3 = gridEdges[(k+1)*ny + j + 1];

    uchar* ec0 = edgeCases + k*ny*(nx-1) + j*(nx-1);
    uchar* ec1 = edgeCases + k*ny*(nx-1) + (j+1)*(nx-1);
    uchar* ec2 = edgeCases + (k+1)*ny*(nx-1) + j*(nx-1);
    uchar* ec3 = edgeCases + (k+1)*ny*(nx-1) + (j+1)*(nx-1);

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
        uchar caseId = calcCubeCase(ec0[i], ec1[i], ec2[i], ec3[i]);

        curCubeCases[i] = caseId;

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

    if(!validKernelSize(gridDim, blockDim))
        std::cout << "GAHHHHHHHHHHHH GP2 ENGINE " << __LINE__ << std::endl; // TODO

    pass2gpu_cubeCases<<<gridDim, blockDim>>>(
        nx, ny, nz,
        edgeCases,
        gridEdges,   // modified
        triCounter,  // modified
        cubeCases);  // modified

    // POSSIBLE to do this here TODO
    // cudaFree(edgeCases);

    if(DEBUG)
    {
       std::cout << "MEOWWWWWW " << cudaGetErrorString(cudaGetLastError()) << std::endl;
    }

    if(DEBUG)
    {
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
    }

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
std::cout << "MEOWWWWWW " << cudaGetErrorString(cudaGetLastError()) << std::endl;

    if(DEBUG)
    {
        size_t sz_ge = nx*ny*sizeof(gridEdge);
        gridEdge* hostges = (gridEdge*)malloc(sz_ge);
        auto w = cudaMemcpy(hostges, gridEdges,
                   sz_ge, cudaMemcpyDeviceToHost);
        if(w != cudaSuccess)
        {
            std::cout << "GHASDCFAKSCLKASCKAS:CKASL:CKAS:DLCKASD:" << std::endl;
            std::cout << cudaGetErrorString(w) << std::endl;
        }

        int sumxstart = 0;
        for(int idx = 0; idx != nx*ny; ++idx)
        {
            sumxstart += hostges[idx].xstart;
        }
        std::cout << "sumxstart " << sumxstart << std::endl;
        free(hostges);
    }

}

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
    int k = blockIdx.y * blockDim.y + threadIdx.y;

    // step 1: accumulate individual y thread
    // step 2: calc block sum
    // step 3: __syncthreads
    // step 4: add to individual y thread

    __shared__ int accum[4*FE_BLOCK_WIDTH];

    if(k < nz)
    {
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

        accum[4*threadIdx.y + 0] = accumX;
        accum[4*threadIdx.y + 1] = accumY;
        accum[4*threadIdx.y + 2] = accumZ;
        accum[4*threadIdx.y + 3] = accumTri;
    }

    __syncthreads();

    if(k < nz)
    {
        if(threadIdx.y == 0) // agh!
        {
            for(int idx = 1; idx != blockDim.y; ++idx)
            {
                accum[4*idx + 0] += accum[4*(idx-1) + 0];
                accum[4*idx + 1] += accum[4*(idx-1) + 1];
                accum[4*idx + 2] += accum[4*(idx-1) + 2];
                accum[4*idx + 3] += accum[4*(idx-1) + 3];
            }

            // answer for global accumulation
            blockAccum[4*blockIdx.y + 0] = accum[4*(blockDim.y-1) + 0];
            blockAccum[4*blockIdx.y + 1] =  accum[4*(blockDim.y-1) + 1];
            blockAccum[4*blockIdx.y + 2] =  accum[4*(blockDim.y-1) + 2];
            blockAccum[4*blockIdx.y + 3] = accum[4*(blockDim.y-1) + 3];
        }
    }
    __syncthreads();

    if(threadIdx.y == 0 || k >= nz)
        return;

    bool isEndK = k == nz-1;
    for(int j = 0; j != ny-1; ++j)
    {
        FlyingEdgesAlgorithm::gridEdge& ge = gridEdges[k*ny + j];

        ge.xstart += accum[4*(threadIdx.y-1) + 0];
        ge.ystart += accum[4*(threadIdx.y-1) + 1];
        ge.zstart += accum[4*(threadIdx.y-1) + 2];

        // put z stuff here..
        if(!isEndK)
            triCounter[k*(ny-1) + j] = accum[4*(threadIdx.y-1) + 3];
    }

    FlyingEdgesAlgorithm::gridEdge& ge = gridEdges[k*ny + (ny-1)];
    ge.xstart += accum[4*(threadIdx.y-1) + 0];
    ge.ystart += accum[4*(threadIdx.y-1) + 1];
    ge.zstart += accum[4*(threadIdx.y-1) + 2];
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

    uint3 gridDim = make_uint3(1, numBlocks, 1);
    uint3 blockDim = make_uint3(1, tz, 1);

    int* hostBlockAccum = (int*)malloc(sizeBlocks);
    for(int idx = 0; idx != 4*numBlocks; ++idx)
    {
        hostBlockAccum[idx] = 0;
    }

    int* deviceBlockAccum;
    cudaMalloc(&deviceBlockAccum, sizeBlocks);

    cudaMemcpy(deviceBlockAccum, hostBlockAccum,
                   sizeBlocks, cudaMemcpyHostToDevice);

    // Accumulate values locally

    if(!validKernelSize(gridDim, blockDim))
        std::cout << "GAHHHHHHHHHHHH GP2 ENGINE " << __LINE__ << std::endl; // TODO


    pass3gpu_blockAccum<<<gridDim, blockDim>>>(
        nx, ny, nz,
        triCounter,
        gridEdges,
        deviceBlockAccum);

    cudaMemcpy(hostBlockAccum, deviceBlockAccum,
               sizeBlocks, cudaMemcpyDeviceToHost);

    if(DEBUG)
    {
        std::cout << "ACCUM ";
        for(int idx = 0; idx != 4*numBlocks; ++idx)
        {
            std::cout << hostBlockAccum[idx] << " ";
        }
        std::cout << std::endl;

        cudaDeviceSynchronize();

        std::cout << "MEOWWWWWW " << cudaGetErrorString(cudaGetLastError()) << std::endl;
    }

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

        // TODO
        if(!validKernelSize(gridDim, blockDim))
            std::cout << "GAHHHHHHHHHHHH GP2 ENGINE " << __LINE__ << std::endl; // TODO

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
    numPoints = hostBlockAccum[4*(numBlocks-1) + 0] +
                hostBlockAccum[4*(numBlocks-1) + 1] +
                hostBlockAccum[4*(numBlocks-1) + 2];
    numTris   = hostBlockAccum[4*(numBlocks-1) + 3];

//    cudaMalloc(&points,  3*sizeof(scalar_t)*numPoints);
//    cudaMalloc(&normals, 3*sizeof(scalar_t)*numPoints);
//    cudaMalloc(&tris, 3*sizeof(int)*numTris);

    if(DEBUG)
    {
        std::cout << "numpoints" << numPoints << std::endl;
        std::cout << "numtris"   << numTris   << std::endl;
    }

    // free memory used in this function
    free(hostBlockAccum);
    cudaFree(deviceBlockAccum);

    cudaDeviceSynchronize();

    if(DEBUG)
    {
        std::cout << "MEOWWWWWW " << cudaGetErrorString(cudaGetLastError()) << std::endl;
    }
}

///////////////////////////////////////////////////////////////////////////////
// Pass 4 of the algorithm
///////////////////////////////////////////////////////////////////////////////
__device__
void computeGradient(
    int const& i, int const& j, int const& k,
    int const& nx, int const& ny, int const& nz,
    scalar_t* data,
    scalar_t* spacing,
    scalar_t* point)
{
    scalar_t x0[2];
    scalar_t x1[2];
    scalar_t x2[2];
    scalar_t run[3];

    size_t dataIdx = k*nx*ny + j*nx + i;

    if (i == 0)
    {
        x0[0] = data[dataIdx + 1];
        x0[1] = data[dataIdx];
        run[0] = spacing[0];
    }
    else if (i == (nx - 1))
    {
        x0[0] = data[dataIdx];
        x0[1] = data[dataIdx - 1];
        run[0] = spacing[0];
    }
    else
    {
        x0[0] = data[dataIdx + 1];
        x0[1] = data[dataIdx - 1];
        run[0] = 2 * spacing[0];
    }

    if (j == 0)
    {
        x1[0] = data[dataIdx + nx];
        x1[1] = data[dataIdx];
        run[1] = spacing[1];
    }
    else if (j == (ny - 1))
    {
        x1[0] = data[dataIdx];
        x1[1] = data[dataIdx - nx];
        run[1] = spacing[1];
    }
    else
    {
        x1[0] = data[dataIdx + nx];
        x1[1] = data[dataIdx - ny];
        run[1] = 2 * spacing[1];
    }

    if (k == 0)
    {
        x2[0] = data[dataIdx + nx*ny];
        x2[1] = data[dataIdx];
        run[2] = spacing[2];
    }
    else if (k == (nz - 1))
    {
        x2[0] = data[dataIdx];
        x2[1] = data[dataIdx - nx*ny];
        run[2] = spacing[2];
    }
    else
    {
        x2[0] = data[dataIdx + nx*ny];
        x2[1] = data[dataIdx - nx*ny];
        run[2] = 2 * spacing[2];
    }

    point[0] = (x0[1] - x0[0]) / run[0];
    point[1] = (x1[1] - x1[0]) / run[1];
    point[2] = (x2[1] - x2[0]) / run[2];
}


__device__
void getCubeInfo(
    int i, int j, int k,
    int nx, int ny, int nz,
    scalar_t* pointValues, scalar_t* zeroPos, scalar_t* spacing,
    scalar_t* pointCube, scalar_t* isovalCube, scalar_t* gradCube)
{
    isovalCube[0] = pointValues[k*ny*nx + j*nx + i];
    isovalCube[1] = pointValues[k*ny*nx + j*nx + i+1];
    isovalCube[2] = pointValues[k*ny*nx + (j+1)*nx + i+1];
    isovalCube[3] = pointValues[k*ny*nx + (j+1)*nx + i];
    isovalCube[4] = pointValues[(k+1)*ny*nx + j*nx + i];
    isovalCube[5] = pointValues[(k+1)*ny*nx + j*nx + i+1];
    isovalCube[6] = pointValues[(k+1)*ny*nx + (j+1)*nx + (i+1)];
    isovalCube[7] = pointValues[(k+1)*ny*nx + (j+1)*nx + i];

    scalar_t xpos = zeroPos[0] + i * spacing[0];
    scalar_t ypos = zeroPos[1] + j * spacing[1];
    scalar_t zpos = zeroPos[2] + k * spacing[2];

    pointCube[0*3 + 0] = xpos;
    pointCube[0*3 + 1] = ypos;
    pointCube[0*3 + 2] = zpos;

    pointCube[1*3 + 0] = xpos + spacing[0];
    pointCube[1*3 + 1] = ypos;
    pointCube[1*3 + 2] = zpos;

    pointCube[2*3 + 0] = xpos + spacing[0];
    pointCube[2*3 + 1] = ypos + spacing[1];
    pointCube[2*3 + 2] = zpos;

    pointCube[3*3 + 0] = xpos;
    pointCube[3*3 + 1] = ypos + spacing[1];
    pointCube[3*3 + 2] = zpos;

    pointCube[4*3 + 0] = xpos;
    pointCube[4*3 + 1] = ypos;
    pointCube[4*3 + 2] = zpos + spacing[2];

    pointCube[5*3 + 0] = xpos + spacing[0];
    pointCube[5*3 + 1] = ypos;
    pointCube[5*3 + 2] = zpos + spacing[2];

    pointCube[6*3 + 0] = xpos + spacing[0];
    pointCube[6*3 + 1] = ypos + spacing[1];
    pointCube[6*3 + 2] = zpos + spacing[2];

    pointCube[7*3 + 0] = xpos;
    pointCube[7*3 + 1] = ypos + spacing[1];
    pointCube[7*3 + 2] = zpos + spacing[2];

    computeGradient(i  , j  , k  , nx, ny, nz, pointValues, spacing, gradCube + 3*0);
    computeGradient(i+1, j  , k  , nx, ny, nz, pointValues, spacing, gradCube + 3*1);
    computeGradient(i+1, j+1, k  , nx, ny, nz, pointValues, spacing, gradCube + 3*2);
    computeGradient(i  , j+1, k  , nx, ny, nz, pointValues, spacing, gradCube + 3*3);
    computeGradient(i  , j  , k+1, nx, ny, nz, pointValues, spacing, gradCube + 3*4);
    computeGradient(i+1, j  , k+1, nx, ny, nz, pointValues, spacing, gradCube + 3*5);
    computeGradient(i+1, j+1, k+1, nx, ny, nz, pointValues, spacing, gradCube + 3*6);
    computeGradient(i  , j+1, k+1, nx, ny, nz, pointValues, spacing, gradCube + 3*7);
}

__device__
void interpolate(
    scalar_t const& weight,
    scalar_t* a,
    scalar_t* b,
    scalar_t* out)
{
    out[0] = a[0] + (weight * (b[0] - a[0]));
    out[1] = a[1] + (weight * (b[1] - a[1]));
    out[2] = a[2] + (weight * (b[2] - a[2]));
}

__device__
void interpolateOnCube(
    uchar const& edge,
    scalar_t const& isoval,
    scalar_t* pts,
    scalar_t* isovals,
    scalar_t* out)
{
    uchar i0 = cuda_util::edgeVertices[edge][0];
    uchar i1 = cuda_util::edgeVertices[edge][1];

    scalar_t weight = (isoval - isovals[i0]) / (isovals[i1] - isovals[i0]);
    interpolate(weight, pts + 3*i0, pts + 3*i1, out);
}

__global__
void pass4gpu_pointsAndNormals(
    int nx, int ny, int nz,
    scalar_t* pointValues, scalar_t* zeroPos, scalar_t* spacing,
    scalar_t isoval,
    FlyingEdgesAlgorithm::gridEdge* gridEdges,
    int* triCounter,
    uchar* cubeCases,
    scalar_t* points, scalar_t* normals, int* tris)
{
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    int k = blockIdx.y * blockDim.y + threadIdx.y;

    if(DEBUG)
    {
        if(j == 0 && k == 0)
        {
//            for(int i = 0; i != 3*1370424; ++i)
//            {
//                points[i] = -1;
//                normals[i] = -1;
//            }

            for(int i = 0; i != 3*2740864; ++i)
                tris[i] = -1;
        }
    }

    if(j >= ny-1 || k >= nz-1)
        return;

    FlyingEdgesAlgorithm::gridEdge& ge0 = gridEdges[k*ny + j];
    FlyingEdgesAlgorithm::gridEdge& ge1 = gridEdges[k*ny + j+1];
    FlyingEdgesAlgorithm::gridEdge& ge2 = gridEdges[(k+1)*ny + j];
    FlyingEdgesAlgorithm::gridEdge& ge3 = gridEdges[(k+1)*ny + j+1];

    int xl, xr;
    calcTrimValues(xl, xr, ge0, ge1, ge2, ge3);

    if(xl == xr)
        return;

    size_t triIdx = triCounter[k*(ny-1) + j];
    uchar* curCubeCaseIds = cubeCases + (nx-1)*(k*(ny-1) + j);

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

    scalar_t pointCube[8*3];
    scalar_t isovalCube[8];
    scalar_t gradCube[8*3];

    for(size_t i = xl; i != xr; ++i)
    {
        bool isXEnd = (i == nx-2);

        uchar caseId = curCubeCaseIds[i];

        if(caseId == 0 || caseId == 255)
        {
            continue;
        }

        const bool* isCut = cuda_util::isCut[caseId]; // has 12 elements

        // Most of the information contained in pointCube, isovalCube
        // and gradCube will be used--but not necessarily all. It has
        // not been tested whether or not obtaining only the information
        // needed will provide a significant speedup--but
        // most likely not.

        // fill out pointCube, isovalCube and gradCube
        getCubeInfo(i, j, k,
                    nx, ny, nz,
                    pointValues, zeroPos, spacing,
                    pointCube, isovalCube, gradCube);

        // Add Points and normals.
        // Calculate global indices for triangles
        int globalIdxs[12];
        if(isCut[0])
        {
            int idx = ge0.xstart + x0counter;
            interpolateOnCube(0, isoval, pointCube, isovalCube, points + 3*idx);
            interpolateOnCube(0, isoval, gradCube, isovalCube, normals + 3*idx);
            globalIdxs[0] = idx;
            ++x0counter;
        }

        if(isCut[3])
        {
            int idx = ge0.ystart + y0counter;
            interpolateOnCube(3, isoval, pointCube, isovalCube, points + 3*idx);
            interpolateOnCube(3, isoval, gradCube, isovalCube, normals + 3*idx);
            globalIdxs[3] = idx;
            ++y0counter;
        }

        if(isCut[8])
        {
            int idx = ge0.zstart + z0counter;
            interpolateOnCube(8, isoval, pointCube, isovalCube, points + 3*idx);
            interpolateOnCube(8, isoval, gradCube, isovalCube, normals + 3*idx);
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
            int idx = ge0.ystart + y0counter;
            if(isXEnd)
            {
                interpolateOnCube(1, isoval, pointCube, isovalCube, points + 3*idx);
                interpolateOnCube(1, isoval, gradCube, isovalCube, normals + 3*idx);
                // y0counter counter doesn't need to be incremented
                // because it won't be used again.
            }
            globalIdxs[1] = idx;
        }

        if(isCut[9])
        {
            int idx = ge0.zstart + z0counter;
            if(isXEnd)
            {
                interpolateOnCube(9, isoval, pointCube, isovalCube, points + 3*idx);
                interpolateOnCube(9, isoval, gradCube, isovalCube, normals + 3*idx);
                // z0counter doesn't need to in incremented.
            }
            globalIdxs[9] = idx;
        }

        if(isCut[2])
        {
            int idx = ge1.xstart + x1counter;
            if(isYEnd)
            {
                interpolateOnCube(2, isoval, pointCube, isovalCube, points + 3*idx);
                interpolateOnCube(2, isoval, gradCube, isovalCube, normals + 3*idx);
            }
            globalIdxs[2] = idx;
            ++x1counter;
        }

        if(isCut[10])
        {
            int idx = ge1.zstart + z1counter;
            if(isYEnd)
            {
                interpolateOnCube(10, isoval, pointCube, isovalCube, points + 3*idx);
                interpolateOnCube(10, isoval, gradCube, isovalCube, normals + 3*idx);
            }
            globalIdxs[10] = idx;
            ++z1counter;
        }

        if(isCut[4])
        {
            int idx = ge2.xstart + x2counter;
            if(isZEnd)
            {
                interpolateOnCube(4, isoval, pointCube, isovalCube, points + 3*idx);
                interpolateOnCube(4, isoval, gradCube, isovalCube, normals + 3*idx);
            }
            globalIdxs[4] = idx;
            ++x2counter;
        }

        if(isCut[7])
        {
            int idx = ge2.ystart + y2counter;
            if(isZEnd)
            {
                interpolateOnCube(7, isoval, pointCube, isovalCube, points + 3*idx);
                interpolateOnCube(7, isoval, gradCube, isovalCube, normals + 3*idx);
            }
            globalIdxs[7] = idx;
            ++y2counter;
        }

        if(isCut[11])
        {
            int idx = ge1.zstart + z1counter;
            if(isXEnd and isYEnd)
            {
                interpolateOnCube(11, isoval, pointCube, isovalCube, points + 3*idx);
                interpolateOnCube(11, isoval, gradCube, isovalCube, normals + 3*idx);
                // z1counter does not need to be incremented.
            }
            globalIdxs[11] = idx;
        }

        if(isCut[5])
        {
            int idx = ge2.ystart + y2counter;
            if(isXEnd and isZEnd)
            {
                interpolateOnCube(5, isoval, pointCube, isovalCube, points + 3*idx);
                interpolateOnCube(5, isoval, gradCube, isovalCube, normals + 3*idx);
                // y2 counter does not need to be incremented.
            }
            globalIdxs[5] = idx;
        }

        if(isCut[6])
        {
            int idx = ge3.xstart + x3counter;
            if(isYEnd and isZEnd)
            {
                interpolateOnCube(6, isoval, pointCube, isovalCube, points + 3*idx);
                interpolateOnCube(6, isoval, gradCube, isovalCube, normals + 3*idx);
            }
            globalIdxs[6] = idx;
            ++x3counter;
        }

        // Add triangles
        const char* caseTri = cuda_util::caseTriangles[caseId]; // size 16
        for(int idx = 0; caseTri[idx] != -1; idx += 3)
        {
            tris[3*triIdx + 0] = i;
            tris[3*triIdx + 1] = j;
            tris[3*triIdx + 2] = k;

//            tris[3*triIdx + 0] = globalIdxs[caseTri[idx]];
//            tris[3*triIdx + 1] = globalIdxs[caseTri[idx+1]];
//            tris[3*triIdx + 2] = globalIdxs[caseTri[idx+2]];
//            ++triIdx;
        }
    }
}


void FlyingEdgesAlgorithm::pass4()
{
    // pass4 calculates points and normals
    //   1) points and normals

    // 1st kernel:           Calculate the main cube rays
    // 2nd and third kernel:

    int ty = 1;//FE_BLOCK_WIDTH_Y / 2; // divide by 2? TODO figure out this problem..
    int tz = 1;//FE_BLOCK_WIDTH_Z / 2; // gah....
    uint3 gridDim = make_uint3(((ny-1) + ty - 1) / ty, ((nz-1) + tz - 1) / tz, 1);
    uint3 blockDim = make_uint3(ty, tz, 1);

    std::cout << gridDim.x << ", " << gridDim.y << ", " << gridDim.z << std::endl;
    std::cout << blockDim.x << ", " << blockDim.y << ", " << blockDim.z << std::endl;

    if(!validKernelSize(gridDim, blockDim))
        std::cout << "GAHHHHHHHHHHHH GP2 ENGINE " << __LINE__ << std::endl; // TODO

    if(DEBUG)
    {
       cudaDeviceSynchronize();
    }

    pass4gpu_pointsAndNormals<<<gridDim, blockDim>>>(
        nx, ny, nz,                                    // input
        pointValues, zeroPos, spacing,                 // input
        isoval,                                        // input
        gridEdges, triCounter, cubeCases,              // input
        points, normals, tris);                        // output

    if(DEBUG)
    {
       cudaDeviceSynchronize();
       std::cout << "MEOWWWWWW " << cudaGetErrorString(cudaGetLastError()) << std::endl;
    }

    if(DEBUG)
    {
        size_t sz = 3 * numPoints * sizeof(scalar_t);

        scalar_t* hostPts = (scalar_t*)malloc(sz);
        scalar_t* hostNrs = (scalar_t*)malloc(sz);
        int*      hostTrs = (int*)malloc(3*numTris*sizeof(int));

        cudaMemcpy(hostPts, points,  sz, cudaMemcpyDeviceToHost);
        cudaMemcpy(hostNrs, normals, sz, cudaMemcpyDeviceToHost);
        cudaMemcpy(hostTrs, tris, 3*numTris*sizeof(int), cudaMemcpyDeviceToHost);

        scalar_t accumP = 0.0;
        for(int idx = 0; idx != 3 * numPoints; ++idx)
        {
            accumP += hostPts[idx];
            accumP += hostTrs[idx];

            while(accumP >= 1000000)
                accumP -= 1000000;
        }

        int accumT = 0;
        int num0 = -1;
        int num9 = 0;
        int num8 = 0;
        int num7 = 0;

        int numSetPoints = 0;
        for(int idx = 0; idx != 3 * numPoints; ++idx)
        {
            if(hostPts[idx] != -1)
                numSetPoints += 1;
        }

        std::cout << "numSetPoints " << numSetPoints << std::endl;

        for(int idx = 0; idx != 3 * numTris; ++idx)
        {
            if(hostTrs[idx] == 0)
                num0 += 1;

            if(hostTrs[idx] == 9)
                num9 += 1;

            if(hostTrs[idx] == 8)
                num8 += 1;

            if(hostTrs[idx] == 7)
                num7 += 1;

            accumT += hostTrs[idx];

            while(accumT >= 1000000)
                accumT -= 1000000;
        }

        std::cout << "pass 4 hashsum " << accumP << ", " << accumT << std::endl;
        std::cout << "num0 in Tris "   << num0 <<                     std::endl;
        std::cout << "num9 in Tris "   << num9 <<                     std::endl;
        std::cout << "num8 in Tris "   << num8 <<                     std::endl;
        std::cout << "num7 in Tris "   << num7 <<                     std::endl;

        for(int idx = 0; idx != numTris*3; idx += 3)
        {
            if(hostTrs[idx] != -1)
            {
                std::cout << hostTrs[idx+0] << ", "
                          << hostTrs[idx+1] << ", "
                          << hostTrs[idx+2] << std::endl;
            }
        }

        free(hostPts);
        free(hostNrs);
        free(hostTrs);
    }
}


