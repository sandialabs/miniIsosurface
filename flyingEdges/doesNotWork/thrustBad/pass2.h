#ifndef PASS2_FE
#define PASS2_FE

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include <thrust/tuple.h>

#include <thrust/functional.h>

#include "Thrust_Config.h"

namespace p2{

struct iter_transform_set_cube_case
  : public thrust::unary_function<
        int,
        thrust::tuple<
            int,
            int,
            typename thrust::device_vector<uchar>::const_pointer,   // ec0
            typename thrust::device_vector<uchar>::const_pointer,   // ec1
            typename thrust::device_vector<uchar>::const_pointer,   // ec2
            typename thrust::device_vector<uchar>::const_pointer,   // ec3
            typename thrust::device_vector<uchar>::pointer> > // cubeCases
{
private:
    using int_ptr_t =   typename thrust::device_vector<int>::const_pointer;
    using const_ptr_t = typename thrust::device_vector<uchar>::const_pointer;
    using ptr_t =       typename thrust::device_vector<uchar>::pointer;
public:
    iter_transform_set_cube_case(
        t3<const int> n,
        int_ptr_t left,
        int_ptr_t right,
        const_ptr_t edgeCases,
        ptr_t cubeCases)
      : n(n),
        left(left),
        right(right),
        edgeCases(edgeCases),
        cubeCases(cubeCases)
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

    __device__
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
        int,
        int,
        const_ptr_t,
        const_ptr_t,
        const_ptr_t,
        const_ptr_t,
        ptr_t>
    operator()(int const& idx) const
    {
        int k = idx / (n.y - 1);
        int j = idx % (n.y - 1);

        int i0 = k*n.y + j;
        int i1 = k*n.y + j + 1;
        int i2 = (k+1)*n.y + j;
        int i3 = (k+1)*n.y + j + 1;

        int left  = calc_left(i0, i1, i2, i3);
        int right = calc_right(i0, i1, i2, i3);

        if(left > right)
            left = right;

        return thrust::make_tuple(
            left,
            right,
            edgeCases + (n.x-1)*i0,
            edgeCases + (n.x-1)*i1,
            edgeCases + (n.x-1)*i2,
            edgeCases + (n.x-1)*i3,
            cubeCases + (n.x-1)*(k*(n.y-1) + j));
    }

    t3<const int> n;
    int_ptr_t left;
    int_ptr_t right;
    const_ptr_t edgeCases;
    ptr_t cubeCases;

    thrust::minimum<int> min_int;
    thrust::maximum<int> max_int;
};

struct set_cube_case
  : public thrust::unary_function<
        thrust::tuple<
            int,
            int,
            typename thrust::device_vector<uchar>::const_pointer,   // ec0
            typename thrust::device_vector<uchar>::const_pointer,   // ec1
            typename thrust::device_vector<uchar>::const_pointer,   // ec2
            typename thrust::device_vector<uchar>::const_pointer,   // ec3
            typename thrust::device_vector<uchar>::pointer>,        // cubeCases
        void>
{
private:
    using const_ptr_t = typename thrust::device_vector<uchar>::const_pointer;
    using ptr_t =       typename thrust::device_vector<uchar>::pointer;

public:
    __device__
    uchar
    calc_cube_case_id(
        uchar const& e0,
        uchar const& e1,
        uchar const& e2,
        uchar const& e3) const
    {
        // e0 | (_,j,k)
        // e1 | (_,j+1,k)
        // e2 | (_,j,k+1)
        // e3 | (_,j+1,k+1)

        uchar caseId = 0;
        if((e0 == 0) || (e0 == 2)) // 0 | (i,j,k)
            caseId |= 1;
        if((e0 == 0) || (e0 == 1)) // 1 | (i+1,j,k)
            caseId |= 2;
        if((e1 == 0) || (e1 == 1)) // 2 | (i+1,j+1,k)
            caseId |= 4;
        if((e1 == 0) || (e1 == 2)) // 3 | (i,j+1,k)
            caseId |= 8;
        if((e2 == 0) || (e2 == 2)) // 4 | (i,j,k+1)
            caseId |= 16;
        if((e2 == 0) || (e2 == 1)) // 5 | (i+1,j,k+1)
            caseId |= 32;
        if((e3 == 0) || (e3 == 1)) // 6 | (i+1,j+1,k+1)
            caseId |= 64;
        if((e3 == 0) || (e3 == 2)) // 7 | (i,j+1,k+1)
            caseId |= 128;
        return caseId;
    }

    __device__
    void
    operator()(
        thrust::tuple<
            int,
            int,
            const_ptr_t,                       // ec0
            const_ptr_t,                       // ec1
            const_ptr_t,                       // ec2
            const_ptr_t,                       // ec3
            ptr_t> tup) const
    {
        int xl = thrust::get<0>(tup);
        int xr = thrust::get<1>(tup);

        const_ptr_t ec0 = thrust::get<2>(tup);
        const_ptr_t ec1 = thrust::get<3>(tup);
        const_ptr_t ec2 = thrust::get<4>(tup);
        const_ptr_t ec3 = thrust::get<5>(tup);

        ptr_t cubeCases = thrust::get<6>(tup);

        for(int i = xl; i != xr; ++i)
        {
            *(cubeCases + i) = calc_cube_case_id(
                *(ec0 + i),
                *(ec1 + i),
                *(ec2 + i),
                *(ec3 + i));
        }
    }

};

///////////////////////////////////////////////////////////////////////////////

struct iter_transform_grid_change
  : public thrust::unary_function<
        int,
        int>
{
    iter_transform_grid_change(t3<const int> n)
      : n(n)
    {}

    __device__
    int
    operator()(int const& idx) const
    {
        int k = idx / (n.y-1);
        int j = idx % (n.y-1);

        return k*n.y + j;
    }

    t3<const int> n;
};

///////////////////////////////////////////////////////////////////////////////

struct iter_transform_cube_ray
  : public thrust::unary_function<
        int,
        thrust::tuple<
            int,
            int,
            typename thrust::device_vector<uchar>::const_pointer> >
{
private:
    using int_ptr_t =   typename thrust::device_vector<int>::const_pointer;
    using const_ptr_t = typename thrust::device_vector<uchar>::const_pointer;
public:
    iter_transform_cube_ray(
        t3<const int> n,
        int_ptr_t left,
        int_ptr_t right,
        const_ptr_t cubeCases)
      : n(n),
        left(left),
        right(right),
        cubeCases(cubeCases)
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

    __device__
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
        int,
        int,
        const_ptr_t>
    operator()(int const& idx) const
    {
        int k = idx / (n.y - 1);
        int j = idx % (n.y - 1);

        int i0 = k*n.y + j;
        int i1 = k*n.y + j + 1;
        int i2 = (k+1)*n.y + j;
        int i3 = (k+1)*n.y + j + 1;

        int left  = calc_left(i0, i1, i2, i3);
        int right = calc_right(i0, i1, i2, i3);

        if(left > right)
            left = right;

        return thrust::make_tuple(
            left,
            right,
            cubeCases + (n.x-1)*(k*(n.y-1) + j));
    }

    t3<const int> n;
    int_ptr_t left;
    int_ptr_t right;
    const_ptr_t cubeCases;

    thrust::minimum<int> min_int;
    thrust::maximum<int> max_int;

};

struct calc_xyz_tri
  : public thrust::unary_function<
        thrust::tuple<
            int,
            int,
            typename thrust::device_vector<uchar>::const_pointer>,
        thrust::tuple<
            int,   // xstart
            int,   // ystart
            int,   // zstart
            int> > // triCount
{
    calc_xyz_tri(
        t3<const int> n)
      : n(n)
    {}

    __device__ // can't put host here because of cuda_util constants
    thrust::tuple<
        int, // xstart
        int, // ystart
        int, // zstart
        int> // triCount
    operator()(thrust::tuple<
        int,
        int,
        typename thrust::device_vector<uchar>::const_pointer> const& tup) const
    {
        int xl = thrust::get<0>(tup);
        int xr = thrust::get<1>(tup);
        auto cubeCases = thrust::get<2>(tup);

        int xstart = 0;
        int ystart = 0;
        int zstart = 0;
        int triCount = 0;

        for(int i = xl; i != xr; ++i)
        {
            uchar case_id = *(cubeCases + i);

            if(case_id == 0 || case_id == 255)
            {
                continue;
            }

           const bool* isCut = cuda_util::isCut[case_id];

            xstart += isCut[0];
            ystart += isCut[3];
            zstart += isCut[8];

            triCount += cuda_util::numTris[case_id];
        }

        if(xr == n.x-1)
        {
            // It could be the case that *(cubeCaes + n.x-2) is either
            // 0 or 255. This is because the is_cut function checks
            // all x, y and z axis and so it is too generous. The alternate
            // is that the is_cut function doesn't check the y and z axis
            // but then there are cases where every cube is cut on a ray but
            // xl == xr..... So too generous is better.
            const bool* isCut = cuda_util::isCut[*(cubeCases + n.x-2)];
            ystart += isCut[1];
            zstart += isCut[9];
        }

        return thrust::make_tuple(
            xstart,
            ystart,
            zstart,
            triCount);
    }

    t3<const int> n;

};

///////////////////////////////////////////////////////////////////////////////

struct iter_transform_xz_cube_ray
  : public thrust::unary_function<
        int,
        thrust::tuple<
            int,
            int,
            typename thrust::device_vector<uchar>::const_pointer> >
{
private:
    using int_ptr_t =   typename thrust::device_vector<int>::const_pointer;
    using const_ptr_t = typename thrust::device_vector<uchar>::const_pointer;
public:
    iter_transform_xz_cube_ray(
        t3<const int> n,
        int_ptr_t left,
        int_ptr_t right,
        const_ptr_t cubeCases)
      : n(n),
        left(left),
        right(right),
        cubeCases(cubeCases)
    {}

    __device__
    int
    calc_left(
        int const& i0,
        int const& i1) const
    {
        return min_int(
            *(left + i0),
            *(left + i1));
    }

    __device__
    int
    calc_right(
        int const& i0,
        int const& i1) const
    {
        return max_int(
            *(right + i0),
            *(right + i1));
    }

    __device__
    thrust::tuple<
        int,
        int,
        const_ptr_t>
    operator()(int const& idx) const
    {
        int k = idx;

        int i0 = k*n.y + n.y-1;
        int i1 = (k+1)*n.y + n.y-1;

        int left  = calc_left(i0, i1);
        int right = calc_right(i0, i1);

        if(left > right)
            left = right;

        return thrust::make_tuple(
            left,
            right,
            cubeCases + (n.x-1)*(k*(n.y-1) + n.y-2));
    }

    t3<const int> n;
    int_ptr_t left;
    int_ptr_t right;
    const_ptr_t cubeCases;

    thrust::minimum<int> min_int;
    thrust::maximum<int> max_int;
};

struct calc_xz
  : public thrust::unary_function<
        thrust::tuple<
            int,
            int,
            typename thrust::device_vector<uchar>::const_pointer>,
        thrust::tuple<
            int,     // xstart
            int> >   // zstart
{
    calc_xz(
        t3<const int> n)
      : n(n)
    {}

    __device__ // can't put host here because of cuda_util constants
    thrust::tuple<
        int, // xstart
        int> // zstart
    operator()(thrust::tuple<
        int,
        int,
        typename thrust::device_vector<uchar>::const_pointer> tup) const
    {
        int xl = thrust::get<0>(tup);
        int xr = thrust::get<1>(tup);
        auto cubeCases = thrust::get<2>(tup);

        int xstart = 0;
        int zstart = 0;


        for(int i = xl; i != xr; ++i)
        {
            uchar case_id = *(cubeCases + i);

            if(case_id == 0 || case_id == 255)
            {
                continue;
            }

            const bool* isCut = cuda_util::isCut[case_id];

            xstart += isCut[2];
            zstart += isCut[10];
        }

        if(xr == n.x-1)
        {
            const bool* isCut = cuda_util::isCut[*(cubeCases + n.x - 2)];
            zstart += isCut[11];
        }

        return thrust::make_tuple(
            xstart,
            zstart);
    }

    t3<const int> n;
};

///////////////////////////////////////////////////////////////////////////////

struct iter_transform_xy_cube_ray
  : public thrust::unary_function<
        int,
        thrust::tuple<
            int,
            int,
            typename thrust::device_vector<uchar>::const_pointer> >
{
private:
    using int_ptr_t =   typename thrust::device_vector<int>::const_pointer;
    using const_ptr_t = typename thrust::device_vector<uchar>::const_pointer;
public:
    iter_transform_xy_cube_ray(
        t3<const int> n,
        int_ptr_t left,
        int_ptr_t right,
        const_ptr_t cubeCases)
      : n(n),
        left(left),
        right(right),
        cubeCases(cubeCases)
    {}

    __device__
    int
    calc_left(
        int const& i0,
        int const& i1) const
    {
        return min_int(
            *(left + i0),
            *(left + i1));
    }

    __device__
    int
    calc_right(
        int const& i0,
        int const& i1) const
    {
        return max_int(
            *(right + i0),
            *(right + i1));
    }

    __device__
    thrust::tuple<
        int,
        int,
        const_ptr_t>
    operator()(int const& idx) const
    {
        int j = idx;

        int i0 = (n.z-1)*n.y + j;
        int i1 = (n.z-1)*n.y + j + 1;

        int left  = calc_left(i0, i1);
        int right = calc_right(i0, i1);

        if(left > right)
            left = right;

        return thrust::make_tuple(
            left,
            right,
            cubeCases + (n.x-1)*((n.z-2)*(n.y-1) + j));
    }

    t3<const int> n;
    int_ptr_t left;
    int_ptr_t right;
    const_ptr_t cubeCases;

    thrust::minimum<int> min_int;
    thrust::maximum<int> max_int;
};

struct calc_xy
  : public thrust::unary_function<
        thrust::tuple<
            int,
            int,
            typename thrust::device_vector<uchar>::const_pointer>,
        thrust::tuple<
            int,     // xstart
            int> >   // ystart
{
    calc_xy(
        t3<const int> n)
      : n(n)
    {}

    __device__ // can't put host here because of cuda_util constants
    thrust::tuple<
        int, // xstart
        int> // ystart
    operator()(thrust::tuple<
        int,
        int,
        typename thrust::device_vector<uchar>::const_pointer> tup) const
    {
        int xl = thrust::get<0>(tup);
        int xr = thrust::get<1>(tup);
        auto cubeCases = thrust::get<2>(tup);

        int xstart = 0;
        int ystart = 0;

        for(int i = xl; i != xr; ++i)
        {
            uchar case_id = *(cubeCases + i);

            if(case_id == 0 || case_id == 255)
            {
                continue;
            }

            const bool* isCut = cuda_util::isCut[case_id];

            xstart += isCut[4];
            ystart += isCut[7];
        }

        if(xr == n.x-1)
        {
            const bool* isCut = cuda_util::isCut[*(cubeCases + n.x-2)];
            ystart += isCut[5];
        }

        return thrust::make_tuple(
            xstart,
            ystart);
    }

    t3<const int> n;
};

///////////////////////////////////////////////////////////////////////////////

struct last_ray
  : public thrust::unary_function<
        int,
        void>
{
private:
    using int_ptr_t =
        typename thrust::device_vector<int>::pointer;
    using int_const_ptr_t =
        typename thrust::device_vector<int>::const_pointer;
    using const_ptr_t =
        typename thrust::device_vector<uchar>::const_pointer;
public:

    last_ray(
        t3<const int> n,
        int_const_ptr_t left,
        int_const_ptr_t right,
        const_ptr_t cubeCases,
        int_ptr_t xstart)
      : n(n),
        left(left),
        right(right),
        cubeCases(cubeCases),
        xstart(xstart)
    {}

    __device__
    void
    operator()(const int& idx) const
    {
        int xl = *(left  + n.y*n.z - 1);
        int xr = *(right + n.y*n.z - 1);

        if(xl > xr)
            xl = xr;

        auto curCubeCases = cubeCases + (n.x-1)*(n.y*n.z - 1);

        int xcount = 0;

        const bool* isCut;

        for(int i = xl; i != xr; ++i)
        {
            uchar case_id = *(curCubeCases + i);

            if(case_id == 0 || case_id == 255)
            {
                continue;
            }

            isCut = cuda_util::isCut[case_id];

            xcount += isCut[6];
        }

        *(xstart + n.y*n.z - 1) = xcount;
    }

    t3<const int> n;
    int_const_ptr_t left;
    int_const_ptr_t right;
    const_ptr_t cubeCases;
    int_ptr_t xstart;
};

}

#endif
