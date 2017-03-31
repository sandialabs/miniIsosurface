#ifndef PASS1_FE
#define PASS1_FE

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include <thrust/tuple.h>

#include <thrust/functional.h>

#include "Thrust_Config.h"

namespace p1{

struct isGE
  : public thrust::unary_function<
        scalar_t,
        bool>
{
    isGE(scalar_t const& isoval)
      : isoval(isoval)
    {}

    __device__
    bool
    operator()(scalar_t const& val)
    {
        return val >= isoval;
    }

    scalar_t const isoval;
};

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
    operator()(int const& idx)
    {
        int k = idx / ((n.x-1)*n.y);
        int j = (idx / (n.x-1)) % n.y;
        int i = idx % (n.x-1);

        return k*n.x*n.y + j*n.x + i;
    }

    t3<const int> n;
};

struct calc_edge_case
  : public thrust::unary_function<
        thrust::tuple<bool, bool>,
        uchar>
{
    __device__
    uchar operator()(thrust::tuple<bool, bool> const& bls) const
    {
        bool const& prevEdge = thrust::get<0>(bls);
        bool const& currEdge = thrust::get<1>(bls);

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
};

struct calc_left_right
  : public thrust::unary_function<
        int,
        thrust::tuple<
            int,
            int> >
{
    calc_left_right(
        t3<const int> n,
        typename thrust::device_vector<uchar>::const_pointer edgeCases)
      : n(n),
        edgeCases(edgeCases)
    {}

    __device__
    thrust::tuple<
        int,
        int>
    operator()(int const& idx) const
    {
        int k = idx / n.y;
        int j = idx % n.y;

        int left = n.x;
        int right = 0;

        for(int i = 0; i != n.x-1; ++i)
        {
            if(is_cut(i, j, k))
            {
                left = i;
                break;
            }
        }

        if(left != n.x)
        {
            for(int i = n.x-2; i != -1; --i)
            {
                if(is_cut(i, j, k))
                {
                    right = i+1;
                    break;
                }
            }
        }

        return thrust::make_tuple(
            left,
            right);
    }

    __device__
    bool is_cut(
        int const& i,
        int const& j,
        int const& k) const
    {
        int edgeCaseIdx = k*(n.x-1)*n.y + j*(n.x-1) + i;
        uchar const& edgeCaseX = *(edgeCases + edgeCaseIdx);
        if(edgeCaseX == 1 || edgeCaseX == 2)
        {
            return true;
        }

        if(j != n.y - 1)
        {
            int edgeCaseIdxY = k*(n.x-1)*n.y + (j+1)*(n.x-1) + i;
            uchar const& edgeCaseY = *(edgeCases + edgeCaseIdxY);

            // If (edgeCaseX, edgeCaseY) is (0, 1), (1, 2), (2, 3), (0, 3)
            //                              (1, 0), (2, 1), (3, 2), (3, 0)
            // and not the other options of (0, 2), (1, 3),
            //                              (2, 0), (3, 1)
            // then the edge along the y axis is cut.
            // So check to see if edgeCaseX + edgeCaseY is odd.
            if((edgeCaseX + edgeCaseY) % 2 == 1)
            {
                return true;
            }
        }

        if(k != n.z - 1)
        {
            size_t edgeCaseIdxZ = (k+1)*(n.x-1)*n.y + j*(n.x-1) + i;
            uchar const& edgeCaseZ = *(edgeCases + edgeCaseIdxZ);

            // Same as above. If it is odd, then there is a cut except this
            // time along the z axis.
            if((edgeCaseX + edgeCaseZ) % 2 == 1)
            {
                return true;
            }
        }

        return false;
    }

    t3<const int> n;
    typename thrust::device_vector<uchar>::const_pointer edgeCases;
};

/*
struct calc_left_right
  : public thrust::unary_function<
        thrust::tuple<
            typename thrust::device_vector<uchar>::const_pointer,
            typename thrust::device_vector<uchar>::const_pointer>,
        thrust::tuple<
            int,
            int> >
{
    calc_left_right(
        t3<const int> n)
      : n(n)
    {}

    __device__
    thrust::tuple<
        int,
        int>
    operator()(
        thrust::tuple<
            typename thrust::device_vector<uchar>::const_pointer,
            typename thrust::device_vector<uchar>::const_pointer> pair) const
    {
        auto it_beg = thrust::get<0>(pair);
        auto it_end = thrust::get<1>(pair);

        int left_side = n.x;
        int right_side = 0;

        for(auto it = it_beg; it != it_end; ++it)
        {
            if (*it == 1 || *it == 2)
            {
                left_side =  thrust::distance(it_beg, it);
                break;
            }
        }

        for(auto it = it_end - 1; it != it_beg; --it)
        {
            if(*it == 1 || *it == 2)
            {
                right_side = thrust::distance(it_beg, it) + 1;
                break;
            }
        }

        // the second for loop stops before once it reaches it_beg
        // TODO use reverse iterator
        if(right_side == 0 && left_side == 0)
        {
            right_side = 1;
        }

        return thrust::make_tuple(
            left_side,
            right_side);
    }

    t3<const int> n;
};
*/

}

#endif
