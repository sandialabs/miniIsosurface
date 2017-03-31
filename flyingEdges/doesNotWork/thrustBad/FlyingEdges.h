#ifndef FLYINGEDGES_H
#define FLYINGEDGES_H

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include <thrust/for_each.h>
#include <thrust/transform.h>

#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>

#include <thrust/execution_policy.h>

#include "../cuda/CudaMarchingCubesTables.h" // TODO

#include "Thrust_Config.h"

#include "util.h"
#include "pass1.h"
#include "pass2.h"
#include "pass3.h" // Lol not even needed TODO
#include "pass4.h"

struct FlyingEdges
{
    FlyingEdges(Image3D const& image, scalar_t const& isoval)
      : image(image),
        isoval(isoval),
        n(image.n),
        left(n.y*n.z),  // not setting
        right(n.y*n.z), // here lar.
        edgeCases((n.x-1)*n.y*n.z),
        cubeCases((n.x-1)*(n.y-1)*(n.z-1)),
        xstart(n.y*n.z),
        ystart(n.y*n.z),
        zstart(n.y*n.z),
        triCount((n.y-1)*(n.z-1)),
        points(0),
        normals(0),
        tris(0)
    {}

    void pass1()
    {
        ///////////////////////////////////////////////////////////////////////
        // set values of trim.edgeCases
        ///////////////////////////////////////////////////////////////////////
        auto isGE_iter =
            thrust::make_transform_iterator(
                image.begin(),
                p1::isGE(isoval));

/*
        int sz = 262144*100;
        std::cout << sz <<  ", " << n.x*n.y*n.z << std::endl;
        std::cout << n.x << " ! " << n.y << " ! " << n.z << std::endl;
        thrust::device_vector<bool> out(sz);
        thrust::copy(
            thrust::device,
            isGE_iter,
            isGE_iter + sz,
            out.begin());
        thrust::host_vector<bool> outh(out);

        int cccount = 0;
        for(bool const& val: outh)
        {
            if(val)
                cccount += 1;
        }
        std::cout << cccount << std::endl;
        return;
*/

        auto element_iterator =
            thrust::make_zip_iterator(
                thrust::make_tuple(
                    isGE_iter,
                    isGE_iter + 1));

        auto index_iterator =
            thrust::make_transform_iterator(
                thrust::make_counting_iterator(0),
                p1::iter_transform_grid_change(n));

        auto input_iter_begin =
            thrust::make_permutation_iterator(
                element_iterator,
                index_iterator);

        auto input_iter_end =
            input_iter_begin + (n.x-1)*n.y*n.z;

        thrust::transform(
            thrust::device,
            input_iter_begin,
            input_iter_end,
            edgeCases.begin(),
            p1::calc_edge_case());

        ///////////////////////////////////////////////////////////////////////

        int count = 0;
        int count0 = 0;
        int count1 = 0;
        int count2 = 0;
        int count3 = 0;

        thrust::host_vector<uchar> host_edgeCases(edgeCases);

        for(uchar const& c: host_edgeCases)
        {
            if(c == 0)
                count0 += 1;
            if(c == 1)
                count1 += 1;
            if(c == 2)
                count2 += 1;
            if(c == 3)
                count3 += 1;
            count += c;
        }
        std::cout << "sum of edgeCases " << count << std::endl;
        std::cout << "counts " << count0 << ", " << count1 << ", "
                               << count2 << ", " << count3 << std::endl;

        ///////////////////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////////////////
        // set values of xl and xr
        ///////////////////////////////////////////////////////////////////////

        auto lr_iterator =
            thrust::make_zip_iterator(
                thrust::make_tuple(
                    left.begin(),
                    right.begin()));

        thrust::transform(
            thrust::device,
            thrust::make_counting_iterator(0),
            thrust::make_counting_iterator(n.y*n.z),
            lr_iterator,
            p1::calc_left_right(
                n,
                edgeCases.data()));

        ///////////////////////////////////////////////////////////////////////

        int countLeft = 0;
        thrust::host_vector<int> host_left(left);
        for(int const& val: host_left)
        {
            if(val == -1)
                countLeft += n.x;
            else
                countLeft += val;
        }

        int countRight = 0;
        thrust::host_vector<int> host_right(right);
        for(int const& val: host_right)
        {
            if(val == -1)
                countRight += 0;
            else
                countRight += val;
        }

        std::cout
            << "xl, xr: " << countLeft << ", " << countRight
            << std::endl;

        ///////////////////////////////////////////////////////////////////////
        cudaDeviceSynchronize();
    }

    void pass2()
    {
        ///////////////////////////////////////////////////////////////////////
        // calculate cube case Ids
        ///////////////////////////////////////////////////////////////////////

        // (nx-1)*(ny-1)*(nz-1) cubes,
        // given i, j, k
        // get [ec(j, k), ec(j+1, k), ec(j, k+1), ec(j+1, k+1)]
        auto edge_case_iterator_begin =
            thrust::make_transform_iterator(
                thrust::make_counting_iterator(0),
                p2::iter_transform_set_cube_case(
                    n,
                    left.data(),
                    right.data(),
                    edgeCases.data(),
                    cubeCases.data()));

        auto edge_case_iterator_end =
            edge_case_iterator_begin + (n.y-1)*(n.z-1);

        thrust::for_each(
            thrust::device,
            edge_case_iterator_begin,
            edge_case_iterator_end,
            p2::set_cube_case());

        ///////////////////////////////////////////////////////////////////////

        thrust::host_vector<uchar> host_cubeCases(cubeCases);
        int count = 0;
        for(uchar const& val: host_cubeCases)
        {
            if(val != 255)
                count += val;
        }

        std::cout << "Cube cases count " << count << std::endl;


        ///////////////////////////////////////////////////////////////////////
        // set xyz start on cubes,
        // set xz start on x, z face
        // set xy start on x, y face
        // set xyz ray
        ///////////////////////////////////////////////////////////////////////

        // xl
        // xr
        // ray of cubeCases
        auto cube_ray_iterator =
            thrust::make_transform_iterator(
                thrust::make_counting_iterator(0),
                p2::iter_transform_cube_ray(
                    n,
                    left.data(),
                    right.data(),
                    cubeCases.data()));

        auto x_iter_grid_change =
            thrust::make_permutation_iterator(
                xstart.begin(),
                thrust::make_transform_iterator(
                    thrust::make_counting_iterator(0),
                    p2::iter_transform_grid_change(n)));

        auto y_iter_grid_change =
            thrust::make_permutation_iterator(
                ystart.begin(),
                thrust::make_transform_iterator(
                    thrust::make_counting_iterator(0),
                    p2::iter_transform_grid_change(n)));

        auto z_iter_grid_change =
            thrust::make_permutation_iterator(
                zstart.begin(),
                thrust::make_transform_iterator(
                    thrust::make_counting_iterator(0),
                    p2::iter_transform_grid_change(n)));

        auto xyz_tri_iterator =
            thrust::make_zip_iterator(
                thrust::make_tuple(
                    x_iter_grid_change,
                    y_iter_grid_change,
                    z_iter_grid_change,
                    triCount.begin()));

        thrust::transform(
            thrust::device,
            cube_ray_iterator,
            cube_ray_iterator + (n.y-1)*(n.z-1),
            xyz_tri_iterator,
            p2::calc_xyz_tri(n));

        thrust::host_vector<int> hostq_xstart(xstart);
        thrust::host_vector<int> hostq_ystart(ystart);
        thrust::host_vector<int> hostq_zstart(zstart);
        thrust::host_vector<int> host_meow(triCount);

        int countqx = 0;
        int countqy = 0;
        int countqz = 0;
        int count_meow = 0;

        for(int i = 0; i != n.y*n.z; ++i)
        {
            countqx += hostq_xstart[i];
            countqy += hostq_ystart[i];
            countqz += hostq_zstart[i];
        }

        for(int i = 0; i != (n.y-1)*(n.z-1); ++i)
        {
            count_meow += host_meow[i];
        }

        std::cout
            << "xyz start: "
            << countqx << ", "
            << countqy << ", "
            << countqz << ", "
            << count_meow << std::endl;

        ///////////////////////////////////////////////////////////////////////

        auto xz_cube_ray_iterator =
            thrust::make_transform_iterator(
                thrust::make_counting_iterator(0),
                p2::iter_transform_xz_cube_ray(
                    n,
                    left.data(),
                    right.data(),
                    cubeCases.data())); // goes from 0 to n.z-1

        auto xz_element_iterator =
            thrust::make_zip_iterator(
                thrust::make_tuple(
                    xstart.begin() + (0*n.y + n.y-1),
                    zstart.begin() + (0*n.y + n.y-1) ));

        auto xz_index_iterator =
            thrust::make_transform_iterator(
                thrust::make_counting_iterator(0),
                util::step(n.y));

        auto xz_iterator =
            thrust::make_permutation_iterator(
                xz_element_iterator,
                xz_index_iterator);

        thrust::transform(
            thrust::device,
            xz_cube_ray_iterator,
            xz_cube_ray_iterator + n.z-1,
            xz_iterator,
            p2::calc_xz(n));

        ///////////////////////////////////////////////////////////////////////

        auto xy_cube_ray_iterator =
            thrust::make_transform_iterator(
                thrust::make_counting_iterator(0),
                p2::iter_transform_xy_cube_ray(
                    n,
                    left.data(),
                    right.data(),
                    cubeCases.data())); // goes from 0 to n.y-1

        auto xy_iterator =
            thrust::make_zip_iterator(
                thrust::make_tuple(
                    xstart.begin() + ((n.z-1)*n.y + 0),
                    ystart.begin() + ((n.z-1)*n.y + 0)));

        thrust::transform(
            thrust::device,
            xy_cube_ray_iterator,
            xy_cube_ray_iterator + n.y-1,
            xy_iterator,
            p2::calc_xy(n));

        ///////////////////////////////////////////////////////////////////////

        thrust::for_each_n(
            thrust::device,
            thrust::make_counting_iterator(225), // 225 is bogus; just a place holder
            1,
            p2::last_ray(
                n,
                left.data(),
                right.data(),
                cubeCases.data(),
                xstart.data()));

        ///////////////////////////////////////////////////////////////////////

        thrust::host_vector<int> host_xstart(xstart);
        thrust::host_vector<int> host_ystart(ystart);
        thrust::host_vector<int> host_zstart(zstart);

        int countx = 0;
        int county = 0;
        int countz = 0;

        for(int i = 0; i != n.y*n.z; ++i)
        {
            countx += host_xstart[i];
            county += host_ystart[i];
            countz += host_zstart[i];
        }

        std::cout
            << "xyz start: "
            << countx << ", "
            << county << ", "
            << countz << std::endl;

        cudaDeviceSynchronize();
    }

    void pass3()
    {
        int tmp_x = xstart[n.y*n.z - 1]; // These transfer (not very much)
        int tmp_y = ystart[n.y*n.z - 1]; // memory.
        int tmp_z = zstart[n.y*n.z - 1]; //
        int tmp_tri = triCount[(n.y-1)*(n.z-1) - 1];

        thrust::exclusive_scan(
            thrust::device,
            triCount.begin(),
            triCount.end(),
            triCount.begin(),
            0); // 0 offset

        thrust::exclusive_scan(
            thrust::device,
            xstart.begin(),
            xstart.end(),
            xstart.begin(),
            0); // 0 offset

        int numX = tmp_x + xstart[n.y*n.z - 1];

        thrust::exclusive_scan(
            thrust::device,
            ystart.begin(),
            ystart.end(),
            ystart.begin(),
            numX); // use the proper offset

        int numY = tmp_y + ystart[n.y*n.z - 1] - numX;

        thrust::exclusive_scan(
            thrust::device,
            zstart.begin(),
            zstart.end(),
            zstart.begin(),
            numX + numY); // use the proper offset

        int numZ = tmp_z + zstart[n.y*n.z - 1] - numX - numY;

        int numPoints = numX + numY + numZ;
        int numTris = tmp_tri + triCount[(n.y-1)*(n.z-1) - 1];

        // This vector is no longer needed
        edgeCases.resize(0);

        // Reserve space for points, normals and tris
        points.resize(numPoints * 3);
        normals.resize(numPoints * 3);
        tris.resize(numTris * 3);

        std::cout << "numPoints " << numPoints << std::endl;
        std::cout << "numTris   " << numTris   << std::endl;

        cudaDeviceSynchronize();
    }

    void pass4()
    {
        //TODO TODO TODO TODO TODO TODO TODO
        //make iterator contain everything and see if it works
        //
        // iter should contain:
        //   j, k
        //   l, r
        //   image
        //   curCubeCases
        //   <x0,y0,z0,x1,z1,x2,y2 counter, triCount counter>
        //   points,
        //   normals,
        //   tris

        auto input_iterator =
            thrust::make_transform_iterator(
                thrust::make_zip_iterator(
                    thrust::make_tuple(
                        thrust::make_counting_iterator(0),
                        thrust::make_constant_iterator(
                            thrust::make_tuple(
                                left.data(),
                                right.data())),
                        thrust::make_constant_iterator(
                            thrust::make_tuple(
                                n.x,
                                n.y,
                                n.z)),
                        thrust::make_constant_iterator(image.data.data()),
                        thrust::make_constant_iterator(cubeCases.data()),
                        thrust::make_constant_iterator(
                            thrust::make_tuple(
                                xstart.data(),
                                ystart.data(),
                                zstart.data(),
                                triCount.data())),
                        thrust::make_constant_iterator(
                            thrust::make_tuple(
                                points.data(),
                                normals.data(),
                                tris.data())))),
                p4::iter_transform_mega());

        thrust::for_each(
            thrust::device,
            input_iterator,
            input_iterator + (n.y-1)*(n.z-1),
            p4::set_data_mega());

/*

        auto input_iterator =
            thrust::make_transform_iterator(
                thrust::make_counting_iterator(0),
                p4::iter_transform_jk_lr(
                    n,
                    left.data(),
                    right.data()));

        thrust::for_each(
            thrust::device,
            input_iterator,
            input_iterator + (n.z-1)*(n.y-1),
            p4::set_points_normals_tris(
                n,                    // input
                isoval,               // input
                image.spacing,        // input
                image.zeroPos,        // input
                image.data.data(),    // input
                cubeCases.data(),     // input
                xstart.data(),        // input
                ystart.data(),        // input
                zstart.data(),        // input
                triCount.data(),      // input
                points.data(),        // output
                normals.data(),       // output
                tris.data()));        // output
*/
        cudaDeviceSynchronize();
        ///////////////////////////////////////////////////////////////////////

        std::cout << "ASDCSDCASCASCASCASCASDCASCASCASC" << std::endl;

        thrust::host_vector<int> host_tris(tris);
        std::cout << host_tris[0] << ", " << host_tris[1] << ", " << host_tris[2] << std::endl;
        int count7 = 0;
        for(int i = 0; i != host_tris.size(); ++i)
        {
            count7 += host_tris[i];
            count7 = count7 % 500000;
        }

        thrust::host_vector<scalar_t> host_points(points);
        scalar_t mmm = 0.0;
        for(scalar_t const& val: host_points)
        {
            mmm += val;
            while(mmm > 10000)
                mmm -= 10000;
            while(mmm < -10000)
                mmm += 10000;
        }

        std::cout << "EHVHAS " << count7 << ", " << mmm << std::endl;
    }

private:
    Image3D const& image;
    scalar_t const isoval;
    t3<const int> n;

    thrust::device_vector<int> left;
    thrust::device_vector<int> right;
    thrust::device_vector<uchar> edgeCases;
    thrust::device_vector<uchar> cubeCases;

    thrust::device_vector<int> xstart;
    thrust::device_vector<int> ystart;
    thrust::device_vector<int> zstart;
    thrust::device_vector<int> triCount;

    thrust::device_vector<scalar_t> points;
    thrust::device_vector<scalar_t> normals;
    thrust::device_vector<int> tris;
};

#endif
