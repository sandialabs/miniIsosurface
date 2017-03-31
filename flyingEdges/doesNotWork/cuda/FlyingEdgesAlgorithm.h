/*
 * flyingEdgesAlgorithm.h
 *
 *  Created on: Feb 17, 2017
 *      Author: dbourge
 */

#ifndef FLYINGEDGESALGORITHM_H_
#define FLYINGEDGESALGORITHM_H_

#include <vector>
#include <array>

#include "../util/FlyingEdges_Config.h"

#include "../util/Image3D.h"
#include "../util/TriangleMesh.h"

#include <iostream> // TODO

struct FlyingEdgesAlgorithm
{
    FlyingEdgesAlgorithm(util::Image3D const& image, scalar_t const& isoval)
      : isoval(isoval),
        nx(image.xdimension()),
        ny(image.ydimension()),
        nz(image.zdimension()),
        outputAllocated(false),
        deallocated(false)
    {
        cudaError_t cuErr;

        // TODO use cudaMalloc3D and cudaMalloc2D

        cuErr = cudaMalloc(&pointValues, nx*ny*nz*sizeof(scalar_t));
        if(cuErr != cudaSuccess)
        {
            std::cout << "~01" << std::endl;
            throw;
        }

        cuErr = cudaMalloc(&zeroPos, 3*sizeof(scalar_t));
        if(cuErr != cudaSuccess)
        {
            std::cout << "~01.3" << std::endl;
            throw;
        }

        cuErr = cudaMalloc(&spacing, 3*sizeof(scalar_t));
        if(cuErr != cudaSuccess)
        {
            std::cout << "~01.5" << std::endl;
            throw;
        }

        cuErr = cudaMalloc(&gridEdges, ny*nz*sizeof(gridEdge));
        if(cuErr != cudaSuccess)
        {
            std::cout << "~02" << std::endl;
            throw;
        }

        cuErr = cudaMalloc(&triCounter, (ny-1)*(nz-1)*sizeof(int));
        if(cuErr != cudaSuccess)
        {
            std::cout << "~03" << std::endl;
            throw;
        }

        cuErr = cudaMalloc(&edgeCases, (nx-1)*ny*nz*sizeof(uchar));
        if(cuErr != cudaSuccess)
        {
            std::cout << "~04" << std::endl;
            throw;
        }

        cuErr = cudaMalloc(&cubeCases, (nx-1)*(ny-1)*(nz-1)*sizeof(uchar));
        if(cuErr != cudaSuccess)
        {
            std::cout << "~05" << std::endl;
            throw;
        }

        // Move image memory onto device
        cudaMemcpy(
            pointValues,
            image.pointer(),
            nx*ny*nz*sizeof(scalar_t),
            cudaMemcpyHostToDevice);

        scalar_t* tmp = (scalar_t*)malloc(3*sizeof(scalar_t));
        auto imageZp = image.getZeroPos();
        tmp[0] = imageZp[0];
        tmp[1] = imageZp[1];
        tmp[2] = imageZp[2];

        cudaMemcpy(
            zeroPos,
            tmp,
            3*sizeof(scalar_t),
            cudaMemcpyHostToDevice);

        auto imageSp = image.getSpacing();
        tmp[0] = imageSp[0];
        tmp[1] = imageSp[1];
        tmp[2] = imageSp[2];

        cudaMemcpy(
            spacing,
            tmp,
            3*sizeof(scalar_t),
            cudaMemcpyHostToDevice);

    }

    ~FlyingEdgesAlgorithm()
    {
        deallocate();
    }

    void deallocate()
    {
        if(!deallocated)
        {
            cudaFree(pointValues);
            cudaFree(zeroPos);
            cudaFree(spacing);
            cudaFree(gridEdges);
            cudaFree(triCounter);
            cudaFree(edgeCases);
            cudaFree(cubeCases);
            if(outputAllocated)
            {
                cudaFree(points);
                cudaFree(normals);
                cudaFree(tris);
            }
        }
        deallocated = true;
    }

    void pass1();

    void pass2();

    void pass3();

    void pass4();

    util::TriangleMesh moveOutput()
    {
        // TODO
        deallocate();
        return util::TriangleMesh();
    }

public:
    struct gridEdge
    {
        gridEdge()
          : xl(0),
            xr(0),
            xstart(0),
            ystart(0),
            zstart(0)
        {}

        // trim values
        // set on pass 1
        int xl;
        int xr;

        // modified on pass 2
        // set on pass 3
        int xstart;
        int ystart;
        int zstart;
    };

private:
    scalar_t* pointValues; // size of nx*ny*nz, the input
    scalar_t* zeroPos;     // size of 3
    scalar_t* spacing;     // size of 3
    scalar_t const isoval;

    int const nx; //
    int const ny; // for indexing
    int const nz; //

    gridEdge* gridEdges; // size of ny*nz
    int* triCounter;     // size of (ny-1)*(nz-1)

    uchar* edgeCases; // size of (nx-1)*ny*nz
    uchar* cubeCases; // size of (nx-1)*(ny-1)*(nz-1)

    bool outputAllocated;
    bool deallocated;

    scalar_t* points;     //
    scalar_t* normals;    //
    int* tris;            // the output
    int numPoints;        //
    int numTris;          //

private:
/*
    inline uchar
    calcCaseEdge(bool const& prevEdge, bool const& currEdge) const;

    inline uchar
    calcCubeCase(uchar const& ec0, uchar const& ec1,
                 uchar const& ec2, uchar const& ec3) const;

    inline void calcTrimValues(
        size_t& xl, size_t& xr, size_t const& j, size_t const& k) const;

    inline std::array<scalar_t, 3>
    interpolateOnCube(
        cube_t const& pts,
        scalarCube_t const& isovals,
        uchar const& edge) const;

    inline std::array<scalar_t, 3>
    interpolate(
        std::array<scalar_t, 3> const& a,
        std::array<scalar_t, 3> const& b,
        scalar_t const& weight) const;
*/
};


#endif
