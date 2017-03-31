#ifndef PASS4_FE
#define PASS4_FE

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include <thrust/tuple.h>

#include <thrust/functional.h>

#include "Thrust_Config.h"

namespace p4 {

struct iter_transform_jk_lr
  : public thrust::unary_function<
        int,
        thrust::tuple<
            int,    // j
            int,    // k
            int,    // left
            int> >  // right
{
    iter_transform_jk_lr(
        t3<const int> n,
        typename thrust::device_vector<int>::const_pointer left,
        typename thrust::device_vector<int>::const_pointer right)
      : n(n),
        left(left),
        right(right)
    {}

    __device__
    int
    calc_left(
        int const& i0,
        int const& i1,
        int const& i2,
        int const& i3) const
    {
        return min_int(
            *(left + i0),
            min_int(
                *(left + i1),
                min_int(
                    *(left + i2),
                    *(left + i3))));
    }

    __device__ // TODO put in util
    int
    calc_right(
        int const& i0,
        int const& i1,
        int const& i2,
        int const& i3) const
    {
        return max_int(
            *(right + i0),
            max_int(
                *(right + i1),
                max_int(
                    *(right + i2),
                    *(right + i3))));
    }

    __device__
    thrust::tuple<
        int,    // j
        int,    // k
        int,    // left
        int>    // right
    operator()(int const& idx) const
    {
        int k = idx / (n.y-1);
        int j = idx % (n.y-1);

        int i0 = k*n.y + j;
        int i1 = k*n.y + j + 1;
        int i2 = (k+1)*n.y + j;
        int i3 = (k+1)*n.y + j + 1;

        int left  = calc_left(i0, i1, i2, i3);
        int right = calc_right(i0, i1, i2, i3);

        if(left > right)
            left = right;

        return thrust::make_tuple(
            j,
            k,
            left,
            right);
    }

    t3<const int> n;
    typename thrust::device_vector<int>::const_pointer left;
    typename thrust::device_vector<int>::const_pointer right;

    thrust::minimum<int> min_int;
    thrust::maximum<int> max_int;
};


struct set_points_normals_tris
  : public thrust::unary_function<
        thrust::tuple<
            int,
            int,
            int,
            int>,
        void>
{
    set_points_normals_tris(
        t3<const int> n,
        scalar_t const& isoval,
        t3<const scalar_t> const spacing,
        t3<const scalar_t> const zeroPos,
        typename thrust::device_vector<scalar_t>::const_pointer data,
        typename thrust::device_vector<uchar>::const_pointer cubeCases,
        typename thrust::device_vector<int>::const_pointer xstart,
        typename thrust::device_vector<int>::const_pointer ystart,
        typename thrust::device_vector<int>::const_pointer zstart,
        typename thrust::device_vector<int>::const_pointer triCount,
        typename thrust::device_vector<scalar_t>::pointer points,
        typename thrust::device_vector<scalar_t>::pointer normals,
        typename thrust::device_vector<int>::pointer tris)
      : n(n),
        isoval(isoval),
        spacing(spacing),
        zeroPos(zeroPos),
        data(data),
        cubeCases(cubeCases),
        xstart(xstart),
        ystart(ystart),
        zstart(zstart),
        triCount(triCount),
        points(points),
        normals(normals),
        tris(tris)
    {}

    __device__
    void
    getPointCube(
        int const& i,
        int const& j,
        int const& k,
        scalar_t* pos)
    {
        scalar_t xpos = 0.0;//zeroPos.x + i * spacing.x;
        scalar_t ypos = 0.0;//zeroPos.y + j * spacing.y;
        scalar_t zpos = 0.0;//zeroPos.z + k * spacing.z;

        pos[0*3 + 0] = xpos;
        pos[0*3 + 1] = ypos;
        pos[0*3 + 2] = zpos;

        pos[1*3 + 0] = xpos + spacing.x;
        pos[1*3 + 1] = ypos;
        pos[1*3 + 2] = zpos;

        pos[2*3 + 0] = xpos + spacing.x;
        pos[2*3 + 1] = ypos + spacing.y;
        pos[2*3 + 2] = zpos;

        pos[3*3 + 0] = xpos;
        pos[3*3 + 1] = ypos + spacing.y;
        pos[3*3 + 2] = zpos;

        pos[4*3 + 0] = xpos;
        pos[4*3 + 1] = ypos;
        pos[4*3 + 2] = zpos + spacing.z;

        pos[5*3 + 0] = xpos + spacing.x;
        pos[5*3 + 1] = ypos;
        pos[5*3 + 2] = zpos + spacing.z;

        pos[6*3 + 0] = xpos + spacing.x;
        pos[6*3 + 1] = ypos + spacing.y;
        pos[6*3 + 2] = zpos + spacing.z;

        pos[7*3 + 0] = xpos;
        pos[7*3 + 1] = ypos + spacing.y;
        pos[7*3 + 2] = zpos + spacing.z;
    }

    __device__
    int
    idxer(
        int const& i,
        int const& j,
        int const& k) const
    {
        return k*n.y*n.x + j*n.x + i;
    }

    __device__
    void
    getValsCube(
        int const& i,
        int const& j,
        int const& k,
        scalar_t* vals)
    {
        vals[0] = *(data + idxer(i,     j,     k));
        vals[1] = *(data + idxer(i + 1, j,     k));
        vals[2] = *(data + idxer(i + 1, j + 1, k));
        vals[3] = *(data + idxer(i    , j + 1, k));
        vals[4] = *(data + idxer(i,     j,     k + 1));
        vals[5] = *(data + idxer(i + 1, j,     k + 1));
        vals[6] = *(data + idxer(i + 1, j + 1, k + 1));
        vals[7] = *(data + idxer(i,     j + 1, k + 1));
    }

    __device__
    void
    computeGradient(
        int const& idx,
        scalar_t* grad,
        int const& i,
        int const& j,
        int const& k)
    {
        scalar_t x[3][2];
        scalar_t run[3];

        int dataIdx = idxer(i, j, k);

        if (i == 0)
        {
            x[0][0] = *(data + (dataIdx + 1));
            x[0][1] = *(data + dataIdx);
            run[0] = spacing.x;
        }
        else if (i == (n.x - 1))
        {
            x[0][0] = *(data + dataIdx);
            x[0][1] = *(data + (dataIdx - 1));
            run[0] = spacing.x;
        }
        else
        {
            x[0][0] = *(data + (dataIdx + 1));
            x[0][1] = *(data + (dataIdx - 1));
            run[0] = 2 * spacing.x;
        }

        if (j == 0)
        {
            x[1][0] = *(data + (dataIdx + n.x));
            x[1][1] = *(data + (dataIdx));
            run[1] = spacing.y;
        }
        else if (j == (n.y - 1))
        {
            x[1][0] = *(data + dataIdx);
            x[1][1] = *(data + (dataIdx - n.x));
            run[1] = spacing.y;
        }
        else
        {
            x[1][0] = *(data + (dataIdx + n.x));
            x[1][1] = *(data + (dataIdx - n.y));
            run[1] = 2 * spacing.y;
        }

        if (k == 0)
        {
            x[2][0] = *(data + (dataIdx + n.x*n.y));
            x[2][1] = *(data + dataIdx);
            run[2] = spacing.z;
        }
        else if (k == (n.z - 1))
        {
            x[2][0] = *(data + dataIdx);
            x[2][1] = *(data + (dataIdx - n.x*n.y));
            run[2] = spacing.z;
        }
        else
        {
            x[2][0] = *(data + (dataIdx + n.x*n.y));
            x[2][1] = *(data + (dataIdx - n.x*n.y));
            run[2] = 2 * spacing.z;
        }

        grad[3*idx + 0] = (x[0][1] - x[0][0]) / run[0];
        grad[3*idx + 1] = (x[1][1] - x[1][0]) / run[1];
        grad[3*idx + 2] = (x[2][1] - x[2][0]) / run[2];
    }

    __device__
    void
    getGradCube(
        int const& i,
        int const& j,
        int const& k,
        scalar_t* grad)
    {
        computeGradient(0, grad, i,     j,     k);
        computeGradient(1, grad, i + 1, j,     k);
        computeGradient(2, grad, i + 1, j + 1, k);
        computeGradient(3, grad, i    , j + 1, k);
        computeGradient(4, grad, i,     j,     k + 1);
        computeGradient(5, grad, i + 1, j,     k + 1);
        computeGradient(6, grad, i + 1, j + 1, k + 1);
        computeGradient(7, grad, i,     j + 1, k + 1);
    }

    __device__
    scalar_t
    interpolate(
        scalar_t const& a,
        scalar_t const& b,
        scalar_t const& weight) const
    {
        return a + (weight * (b - a));
    }

    __device__
    void interpolateOnCube(
        scalar_t* pts,
        scalar_t* isovals,
        int const& edgeNum,
        typename thrust::device_vector<scalar_t>::pointer out,
        int const& idx)
    {
        uchar i0 = cuda_util::edgeVertices[edgeNum][0];
        uchar i1 = cuda_util::edgeVertices[edgeNum][1];

        scalar_t weight = (isoval - isovals[i0]) / (isovals[i1] - isovals[i0]);

        *(out + 3*idx + 0) = interpolate(pts[3*i0 + 0], pts[3*i1 + 0], weight);
        *(out + 3*idx + 1) = interpolate(pts[3*i0 + 1], pts[3*i1 + 1], weight);
        *(out + 3*idx + 2) = interpolate(pts[3*i0 + 2], pts[3*i1 + 2], weight);
    }

    __device__
    void add_points_and_normals(
        const bool* isCut,
        scalar_t* pointCube,
        scalar_t* isovalCube,
        scalar_t* gradCube,
        int* globalIdxs,
        int const& edgeNum,
        bool const& add,
        int& counter)
    {
        if(isCut[edgeNum])
        {
            if(add)
            {
                interpolateOnCube(
                    pointCube,
                    isovalCube,
                    edgeNum,
                    points,
                    counter);
                interpolateOnCube(
                    gradCube,
                    isovalCube,
                    edgeNum,
                    normals,
                    counter);
            }
            globalIdxs[edgeNum] = counter;
            ++counter;
        }
    }

    __device__
    void
    operator()(
        thrust::tuple<
            int,                            // j
            int,                            // k
            int,                            // left
            int> const& tup)                // right
    {
        int j = thrust::get<0>(tup);
        int k = thrust::get<1>(tup);
        int left  = thrust::get<2>(tup);
        int right = thrust::get<3>(tup);

        int x0counter = *(xstart + (k*n.y + j));
        int y0counter = *(ystart + (k*n.y + j));
        int z0counter = *(zstart + (k*n.y + j));

        int x1counter = *(xstart + (k*n.y + j+1));
        int z1counter = *(zstart + (k*n.y + j+1));

        int x2counter = *(xstart + ((k+1)*n.y + j));
        int y2counter = *(ystart + ((k+1)*n.y + j));

        int x3counter = *(xstart + ((k+1)*n.y + j+1));

        int triIdx = *(tris + (k*(n.y-1) + j));

        bool isYEnd = (j == n.y-2);
        bool isZEnd = (k == n.z-2);

        auto curCubeCases = cubeCases + k*(n.x-1)*(n.y-1) + j*(n.x-1);

        for(int i = left; i != right; ++i)
        {
            bool isXEnd = (i == n.x-2);

            uchar caseId = *(curCubeCases + i);

            if(caseId == 0 || caseId == 255)
            {
                continue;
            }

            const bool* isCut = cuda_util::isCut[caseId];

            scalar_t pointCube[24];
            getPointCube(i, j, k, pointCube);

            scalar_t isovalCube[8];
            getValsCube(i, j, k, isovalCube);

            scalar_t gradCube[24];
            getGradCube(i, j, k, gradCube);

            int globalIdxs[12];

            // add_points_and_normals adds to points and normals
            // if fit; updates globalIdxs and increments counter.
            add_points_and_normals(
                   isCut, pointCube, isovalCube, gradCube, globalIdxs,
                                    0,             true,   x0counter);
            add_points_and_normals(
                   isCut, pointCube, isovalCube, gradCube, globalIdxs,
                                    3,             true,   y0counter);
            add_points_and_normals(
                   isCut, pointCube, isovalCube, gradCube, globalIdxs,
                                    8,             true,   z0counter);

            add_points_and_normals(
                   isCut, pointCube, isovalCube, gradCube, globalIdxs,
                                    1,           isXEnd,   y0counter);
            add_points_and_normals(
                   isCut, pointCube, isovalCube, gradCube, globalIdxs,
                                    9,           isXEnd,   z0counter);

            add_points_and_normals(
                   isCut, pointCube, isovalCube, gradCube, globalIdxs,
                                    2,           isYEnd,   x1counter);
            add_points_and_normals(
                   isCut, pointCube, isovalCube, gradCube, globalIdxs,
                                   10,           isYEnd,   z1counter);

            add_points_and_normals(
                   isCut, pointCube, isovalCube, gradCube, globalIdxs,
                                    4,           isZEnd,   x2counter);
            add_points_and_normals(
                   isCut, pointCube, isovalCube, gradCube, globalIdxs,
                                    7,           isZEnd,   y2counter);

            add_points_and_normals(
                   isCut, pointCube, isovalCube, gradCube, globalIdxs,
                                   11, isXEnd && isYEnd,   z1counter);

            add_points_and_normals(
                   isCut, pointCube, isovalCube, gradCube, globalIdxs,
                                    5, isXEnd && isZEnd,   y2counter);

            add_points_and_normals(
                   isCut, pointCube, isovalCube, gradCube, globalIdxs,
                                    6, isYEnd && isZEnd,   x3counter);

            const char* caseTri = cuda_util::caseTriangles[caseId];
            for(int idx = 0; caseTri[idx] != -1; idx += 3)
            {
                tris[3*triIdx + 0] = globalIdxs[caseTri[idx + 0]];
                tris[3*triIdx + 1] = globalIdxs[caseTri[idx + 1]];
                tris[3*triIdx + 2] = globalIdxs[caseTri[idx + 2]];
                ++triIdx;
            }
        }
    }

    t3<const int> n;
    scalar_t const isoval;
    t3<const scalar_t> spacing;
    t3<const scalar_t> zeroPos;
    typename thrust::device_vector<scalar_t>::const_pointer data;
    typename thrust::device_vector<uchar>::const_pointer cubeCases;
    typename thrust::device_vector<int>::const_pointer xstart;
    typename thrust::device_vector<int>::const_pointer ystart;
    typename thrust::device_vector<int>::const_pointer zstart;
    typename thrust::device_vector<int>::const_pointer triCount;

    typename thrust::device_vector<scalar_t>::pointer points;
    typename thrust::device_vector<scalar_t>::pointer normals;
    typename thrust::device_vector<int>::pointer tris;
};

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//
//
//
//
//
//
//
//
//
//
////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
//
//
///////////////
////
///
///
///
///
///
//
//
//
//
//
//
//
//
//
////////////////////////////
//


template <typename T>
using vec_ptr = typename thrust::device_vector<T>::pointer;
template <typename T>
using const_vec_ptr = typename thrust::device_vector<T>::const_pointer;

struct iter_transform_mega
  : public thrust::unary_function<
        thrust::tuple<
            int,                        // idx
            thrust::tuple<
                const_vec_ptr<int>,
                const_vec_ptr<int> >,
            thrust::tuple<
                int,                    //  nx
                int,                    //  ny
                int>,                   //  nz
            const_vec_ptr<scalar_t>,    // image
            const_vec_ptr<uchar>,       // cubeCases,
            thrust::tuple<
                const_vec_ptr<int>,     // xstart
                const_vec_ptr<int>,     // ystart
                const_vec_ptr<int>,     // zstart
                const_vec_ptr<int> >,   // triCount
            thrust::tuple<
                vec_ptr<scalar_t>,      // points
                vec_ptr<scalar_t>,      // normals
                vec_ptr<int> > >,       // tris
        thrust::tuple<
            thrust::tuple<
                int,           //  j
                int>,          //  k
            thrust::tuple<
                int,           //  left
                int>,          //  right
            thrust::tuple<
                int,           // nx
                int,           // ny
                int>,          // nz
            thrust::tuple<
                int,           // x0     counters...
                int,           // y0
                int,           // z0
                int,           // x1
                int,           // z1
                int,           // x2
                int,           // y2
                int,           // x3
                int>,          // tri
            thrust::tuple<
                const_vec_ptr<scalar_t>,        //image
                const_vec_ptr<uchar> >,         // cubeCases
            thrust::tuple<
                vec_ptr<scalar_t>, // points   output...
                vec_ptr<scalar_t>, // normals
                vec_ptr<int> > > > // tris
{
    __device__
    int
    calc_left(
        const_vec_ptr<int> left,
        int const& i0,
        int const& i1,
        int const& i2,
        int const& i3) const
    {
        return min_int(
            *(left + i0),
            min_int(
                *(left + i1),
                min_int(
                    *(left + i2),
                    *(left + i3))));
    }

    __device__ // TODO put in util
    int
    calc_right(
        const_vec_ptr<int> right,
        int const& i0,
        int const& i1,
        int const& i2,
        int const& i3) const
    {
        return max_int(
            *(right + i0),
            max_int(
                *(right + i1),
                max_int(
                    *(right + i2),
                    *(right + i3))));
    }

    __device__
    thrust::tuple<
        thrust::tuple<
            int,           //  j
            int>,          //  k
        thrust::tuple<
            int,           //  left
            int>,          //  right
        thrust::tuple<
            int,           // nx
            int,           // ny
            int>,          // nz
        thrust::tuple<
            int,           // x0     counters...
            int,           // y0
            int,           // z0
            int,           // x1
            int,           // z1
            int,           // x2
            int,           // y2
            int,           // x3
            int>,          // tri
        thrust::tuple<
            const_vec_ptr<scalar_t>,        //image
            const_vec_ptr<uchar> >,         // cubeCases
        thrust::tuple<
            vec_ptr<scalar_t>, // points   output...
            vec_ptr<scalar_t>, // normals
            vec_ptr<int> > >   // tris
    operator()(
        thrust::tuple<
            int,                        // idx
            thrust::tuple<
                const_vec_ptr<int>,     // leftvec
                const_vec_ptr<int> >,   // rightvec
            thrust::tuple<
                int,                    //  nx
                int,                    //  ny
                int>,                   //  nz
            const_vec_ptr< const scalar_t>,    // image
            const_vec_ptr<uchar>,       // cubeCases,
            thrust::tuple<
                const_vec_ptr<int>,     // xstart
                const_vec_ptr<int>,     // ystart
                const_vec_ptr<int>,     // zstart
                const_vec_ptr<int> >,   // triCount
            thrust::tuple<
                vec_ptr<scalar_t>,      // points
                vec_ptr<scalar_t>,      // normals
                vec_ptr<int> > >&       // tris
        tup) const
    {
        int& idx = thrust::get<0>(tup);

        int& nx = thrust::get<0>(thrust::get<2>(tup));
        int& ny = thrust::get<1>(thrust::get<2>(tup));
        int& nz = thrust::get<2>(thrust::get<2>(tup));

        int k = idx / (ny-1);
        int j = idx % (ny-1);

        int i0 = k*ny + j;
        int i1 = k*ny + j + 1;
        int i2 = (k+1)*ny + j;
        int i3 = (k+1)*ny + j + 1;

        const_vec_ptr<int>& leftvec  = thrust::get<0>(thrust::get<1>(tup));
        const_vec_ptr<int>& rightvec = thrust::get<1>(thrust::get<1>(tup));

        int left  = calc_left(leftvec, i0, i1, i2, i3);
        int right = calc_right(rightvec, i0, i1, i2, i3);

        if(left > right)
            left = right;

        const_vec_ptr<int>& xstart = thrust::get<0>(thrust::get<5>(tup));
        const_vec_ptr<int>& ystart = thrust::get<1>(thrust::get<5>(tup));
        const_vec_ptr<int>& zstart = thrust::get<2>(thrust::get<5>(tup));
        const_vec_ptr<int>& tri =    thrust::get<3>(thrust::get<5>(tup));

        int x0 = *(xstart + (k*ny + j));
        int y0 = *(ystart + (k*ny + j));
        int z0 = *(zstart + (k*ny + j));

        int x1 = *(xstart + (k*ny + j+1));
        int z1 = *(zstart + (k*ny + j+1));

        int x2 = *(xstart + ((k+1)*ny + j));
        int y2 = *(ystart + ((k+1)*ny + j));

        int x3 = *(xstart + ((k+1)*ny + j+1));

        int tr = *(tri + (k*(ny-1) + j));

        const_vec_ptr<scalar_t>& image = thrust::get<3>(tup);
        const_vec_ptr<uchar> curCubeCases =
            thrust::get<4>(tup) + k*(ny-1)*(nx-1) + j*(nx-1);

        vec_ptr<scalar_t>& points  = thrust::get<0>(thrust::get<6>(tup));
        vec_ptr<scalar_t>& normals = thrust::get<1>(thrust::get<6>(tup));
        vec_ptr<int>& triouts      = thrust::get<2>(thrust::get<6>(tup));

        return thrust::make_tuple(
            thrust::make_tuple(
                j,
                k),
            thrust::make_tuple(
                left,
                right),
            thrust::make_tuple(
                nx,
                ny,
                nz),
            thrust::make_tuple(
                x0,
                y0,
                z0,
                x1,
                z1,
                x2,
                y2,
                x3,
                tr),
            thrust::make_tuple(
                image,
                curCubeCases),
            thrust::make_tuple(
                points,
                normals,
                triouts));
    }

    thrust::minimum<int> min_int;
    thrust::maximum<int> max_int;
};



struct set_data_mega
  : public thrust::unary_function<
        thrust::tuple<
            thrust::tuple<
                int,           //  j
                int>,          //  k
            thrust::tuple<
                int,           //  left
                int>,          //  right
            thrust::tuple<
                int,           // nx
                int,           // ny
                int>,          // nz
            thrust::tuple<
                int,           // x0     counters...
                int,           // y0
                int,           // z0
                int,           // x1
                int,           // z1
                int,           // x2
                int,           // y2
                int,           // x3
                int>,          // tr
            thrust::tuple<
                const_vec_ptr<scalar_t>,        // image
                const_vec_ptr<uchar> >,         // curCubeCases
            thrust::tuple<
                vec_ptr<scalar_t>, // points   output...
                vec_ptr<scalar_t>, // normals
                vec_ptr<int> > >, // tris
        void>
{
    __device__
    void
    operator()(
        thrust::tuple<
            thrust::tuple<
                int,           //  j
                int>,          //  k
            thrust::tuple<
                int,           //  left
                int>,          //  right
            thrust::tuple<
                int,           // nx
                int,           // ny
                int>,          // nz
            thrust::tuple<
                int,           // x0     counters...
                int,           // y0
                int,           // z0
                int,           // x1
                int,           // z1
                int,           // x2
                int,           // y2
                int,           // x3
                int>,          // tri
            thrust::tuple<
                const_vec_ptr<scalar_t>,        //image
                const_vec_ptr<uchar> >,     // cubeCases
            thrust::tuple<
                vec_ptr<scalar_t>, // points   output...
                vec_ptr<scalar_t>, // normals
                vec_ptr<int> > >& // tris
         tup) const
    {
        int const& j = thrust::get<0>(thrust::get<0>(tup));
        int const& k = thrust::get<1>(thrust::get<0>(tup));

        int const& left  = thrust::get<0>(thrust::get<1>(tup));
        int const& right = thrust::get<1>(thrust::get<1>(tup));

        int const& nx = thrust::get<0>(thrust::get<2>(tup));
        int const& ny = thrust::get<1>(thrust::get<2>(tup));
        int const& nz = thrust::get<2>(thrust::get<2>(tup));

        int const& x0 = thrust::get<0>(thrust::get<3>(tup));
        int const& y0 = thrust::get<1>(thrust::get<3>(tup));
        int const& z0 = thrust::get<2>(thrust::get<3>(tup));

        int const& x1 = thrust::get<3>(thrust::get<3>(tup));
        int const& z1 = thrust::get<4>(thrust::get<3>(tup));

        int const& x2 = thrust::get<5>(thrust::get<3>(tup));
        int const& y2 = thrust::get<6>(thrust::get<3>(tup));

        int const& x3 = thrust::get<7>(thrust::get<3>(tup));

        int const& xr = thrust::get<8>(thrust::get<3>(tup));

        const_vec_ptr<scalar_t>& image =
            thrust::get<0>(thrust::get<4>(tup));
        const_vec_ptr<uchar>& curCubeCases =
            thrust::get<1>(thrust::get<4>(tup));

        vec_ptr<scalar_t>& points  = thrust::get<0>(thrust::get<5>(tup));
        vec_ptr<scalar_t>& normals = thrust::get<1>(thrust::get<5>(tup));
        vec_ptr<int>& tris         = thrust::get<2>(thrust::get<5>(tup));
    }
};









}

#endif
