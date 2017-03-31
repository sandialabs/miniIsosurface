#include <iostream>
#include <string.h>
#include <cstdlib>

#include "LoadImage.h"
#include "SaveTriangleMesh.h"

#include "../util/Timer.h"
#include "../mantevoCommon/YAML_Doc.hpp"

#include "config.h"
#include "structs.h"

// TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
//
// While correct, this code is slower than the serial version.
// The most likely culprit is that way too much memory is being
// allocated upfront.
//
// Instead of setting points and normals in the first part--and
// thereby allocating way too much memory that is mostly filled with
// unused values, the next version should only fill out a0, b0, c0,
// triscan and caseid on the first part. Then scan them. After the
// scan step, allocate however much memory available and set points
// and normals and then triangles.
//
// The limiting factor is not how much computation can be done but
// memory.
//
// TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO

int main(int argc, char* argv[])
{
    float isoval;
    bool isovalSet = false;
    char* vtkFile = NULL;
    char* outFile = NULL;
    std::string yamlDirectory = "";
    std::string yamlFileName  = "";
    int process_size = 0; // number of elements to process. If not set,
                          // processes one y x z plane at a time.

    // Read command line arguments
    for(int i=0; i<argc; i++)
    {
        if( (strcmp(argv[i], "-i") == 0) || (strcmp(argv[i], "-input_file") == 0))
        {
            vtkFile = argv[++i];
        }
        else if( (strcmp(argv[i], "-o") == 0) || (strcmp(argv[i], "-output_file") == 0))
        {
            outFile = argv[++i];
        }
        else if( (strcmp(argv[i], "-v") == 0) || (strcmp(argv[i], "-isoval") == 0))
        {
            isovalSet = true;
            isoval = atof(argv[++i]);
        }
        else if( (strcmp(argv[i], "-p") == 0) || (strcmp(argv[i], "-process_size") == 0))
        {
            process_size = atoi(argv[++i]);
        }
        else if( (strcmp(argv[i], "-y") == 0) || (strcmp(argv[i], "-yaml_output_file") == 0))
        {
            std::string wholeFile(argv[++i]);

            std::size_t pos = wholeFile.rfind("/");
            if(pos == std::string::npos)
            {
                yamlDirectory = "./";
                yamlFileName = wholeFile;
            }
            else
            {
                yamlDirectory = wholeFile.substr(0, pos + 1);
                yamlFileName = wholeFile.substr(pos + 1);
            }
        }
        else if( (strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "-help") == 0))
        {
            std::cout <<
                "Serial Flying Edges Options:"    << std::endl <<
                "  -input_file (-i)"              << std::endl <<
                "  -output_file (-o)"             << std::endl <<
                "  -isoval (-v)"                  << std::endl <<
                "  -process_size (-p)"            << std::endl <<
                "  -yaml_output_file (-y)"        << std::endl <<
                "  -help (-h)"                    << std::endl;
            return 0;
        }
    }

    if(isovalSet == false || vtkFile == NULL || outFile == NULL)
    {
        std::cout << "Error: isoval, input_file and output_file must be set." << std::endl <<
                     "Try -help" << std::endl;
        return 0;
    }

    std::vector<scalar_t> data_host;
    scalar_t spacing_x, spacing_y, spacing_z;
    scalar_t zeropos_x, zeropos_y, zeropos_z;
    int nx, ny, nz;

    util::loadImage_thrust(
        vtkFile,
        data_host,
        spacing_x, spacing_y, spacing_z,
        zeropos_x, zeropos_y, zeropos_z,
        nx, ny, nz);

    YAML_Doc doc("Modified Flying Edges", "0.1", yamlDirectory, yamlFileName);

    doc.add("Modified Flying Edges Algorithm", "cuda, thrust");
    doc.add("Volume image data file path", vtkFile);
    doc.add("Polygonal mesh output file", outFile);
    doc.add("Isoval", isoval);

    if(process_size == 0)
    {
        process_size = ny*nz;
    }
    doc.add("Process Size", process_size);

    doc.add("File x-dimension", nx);
    doc.add("File y-dimension", ny);
    doc.add("File z-dimension", nz);

    ///////////////////////////////////////////////////////////////////////////
    // copy over image data to device and Allocate memory

    // Start timer
    util::Timer run_time;

    util::Timer run_time_allocate_memory;

    vector<scalar_t> image_data(data_host);

    int n = nx*ny*nz;
    int processed = 0;

    // a_xyz
    vector<scalar_t> ax(process_size);  // points    ; will be set to -1 if not needed
    vector<scalar_t> ay(process_size);  //           ; the others don't matter
    vector<scalar_t> az(process_size);
    vector<scalar_t> axn(process_size); // normals
    vector<scalar_t> ayn(process_size);
    vector<scalar_t> azn(process_size);

    // b_xyz
    vector<scalar_t> bx(process_size); // has to be set to -1 if not used
    vector<scalar_t> by(process_size);
    vector<scalar_t> bz(process_size);
    vector<scalar_t> bxn(process_size);
    vector<scalar_t> byn(process_size);
    vector<scalar_t> bzn(process_size);

    // c_xyz
    vector<scalar_t> cx(process_size); // has to be set to -1 if not used
    vector<scalar_t> cy(process_size);
    vector<scalar_t> cz(process_size);
    vector<scalar_t> cxn(process_size);
    vector<scalar_t> cyn(process_size);
    vector<scalar_t> czn(process_size);

    // cube_ids
    vector<uchar> cube_ids(n);

    // tri_scan
    vector<int> tri_scan(n); // needs to be int for scan step

    vector<int> a0(n);
    vector<int> b0(n);
    vector<int> c0(n);

    int num_triangles = 0;
    int num_points = 0;

    scalar_t default_value = std::min(
        zeropos_x,
        std::min(zeropos_y, zeropos_z)) - 1;

    int max_cur_points = process_size;
    vector<scalar_t> pts_x(process_size); // guess at the size, may make it bigger
    vector<scalar_t> pts_y(process_size); // during algorithm. But probably already
    vector<scalar_t> pts_z(process_size); // too big of a guess

    vector<scalar_t> nrs_x(process_size);
    vector<scalar_t> nrs_y(process_size);
    vector<scalar_t> nrs_z(process_size);

    host_vector<scalar_t> host_pts_x(0);
    host_vector<scalar_t> host_pts_y(0);
    host_vector<scalar_t> host_pts_z(0);
    host_vector<scalar_t> host_nrs_x(0);
    host_vector<scalar_t> host_nrs_y(0);
    host_vector<scalar_t> host_nrs_z(0);

    run_time_allocate_memory.stop();

    util::Timer run_time_points_and_normals;

    while(processed != n)
    {
        int p = std::min(process_size, n - processed);

        auto a_iter = make_zip_iterator(
            make_tuple(
                ax.begin(),  ay.begin(),  az.begin(),
                axn.begin(), ayn.begin(), azn.begin()));
        auto b_iter = make_zip_iterator(
            make_tuple(
                bx.begin(),  by.begin(),  bz.begin(),
                bxn.begin(), byn.begin(), bzn.begin()));
        auto c_iter = make_zip_iterator(
            make_tuple(
                cx.begin(),  cy.begin(),  cz.begin(),
                cxn.begin(), cyn.begin(), czn.begin()));

        auto pts_nors_plus = make_zip_iterator(
            make_tuple(
                a_iter,
                b_iter,
                c_iter,
                cube_ids.begin() + processed,
                tri_scan.begin() + processed));

        transform(
            policy,
            make_counting_iterator(processed),     // Will calculate v0, ..., v7 from
            make_counting_iterator(processed + p), // function. As well as gradient vs
            pts_nors_plus,
            abc_transform(
                nx, ny, nz,
                spacing_x, spacing_y, spacing_z,
                zeropos_x, zeropos_y, zeropos_z,
                isoval, image_data.data()));

        ///////////////////////////////////////////////////////////////////////////

        // tmp sums
        int np_a_temp = (ax[p-1] != -1);
        int np_b_temp = (bx[p-1] != -1);
        int np_c_temp = (cx[p-1] != -1);
        int num_triangles_temp = tri_scan[processed + p - 1];

        auto ax_mod_iter = make_transform_iterator(
            ax.begin(),
            neq(default_value));
        auto bx_mod_iter = make_transform_iterator(
            bx.begin(),
            neq(default_value));
        auto cx_mod_iter = make_transform_iterator(
            cx.begin(),
            neq(default_value));

        exclusive_scan(
            policy,
            ax_mod_iter,
            ax_mod_iter + p,
            a0.begin() + processed,
            num_points);

        exclusive_scan(
            policy,
            bx_mod_iter,
            bx_mod_iter + p,
            b0.begin() + processed,
            np_a_temp + a0[processed + p - 1]); // increase starting value

        exclusive_scan(
            policy,
            cx_mod_iter,
            cx_mod_iter + p,
            c0.begin() + processed,
            np_b_temp + b0[processed + p - 1]); // increase starting value

        // Don't need num triangles; can do based off of cubeIds!
        exclusive_scan(
            policy,
            tri_scan.begin() + processed,
            tri_scan.begin() + processed + p,
            tri_scan.begin() + processed,
            num_triangles);

        int prev_num_points = num_points;
        // final values of sums
        num_triangles = num_triangles_temp + tri_scan[processed + p - 1];
        num_points = np_c_temp + c0[processed + p - 1];

        ///////////////////////////////////////////////////////////////////////////
        // Allocate points and normals. Triangles will be allocated later
        int cur_num_points = num_points - prev_num_points;

        if(cur_num_points > max_cur_points)
        {
            max_cur_points = cur_num_points;

            pts_x.resize(cur_num_points);
            pts_y.resize(cur_num_points);
            pts_z.resize(cur_num_points);

            nrs_x.resize(cur_num_points);
            nrs_y.resize(cur_num_points);
            nrs_z.resize(cur_num_points);
        }

        ///////////////////////////////////////////////////////////////////////////
        // Set points and normals
        auto indexer = make_zip_iterator(
            make_tuple(
                a0.begin() + processed,
                b0.begin() + processed,
                c0.begin() + processed));

        auto indexer_plus_info_iterator = make_zip_iterator(
            make_tuple(
                indexer,
                a_iter,
                b_iter,
                c_iter));

        for_each(
            policy,
            indexer_plus_info_iterator,
            indexer_plus_info_iterator + p,
            set_points_and_normals(
                pts_x.data(), pts_y.data(), pts_z.data(),
                nrs_x.data(), nrs_y.data(), nrs_z.data(),
                default_value,
                prev_num_points));

        host_pts_x.resize(num_points);
        host_pts_y.resize(num_points);
        host_pts_z.resize(num_points);
        host_nrs_x.resize(num_points);
        host_nrs_y.resize(num_points);
        host_nrs_z.resize(num_points);

        int& cnp = cur_num_points; // just making the names shorter
        int& pnp = prev_num_points;

        thrust::copy(pts_x.begin(), pts_x.begin() + cnp, host_pts_x.begin() + pnp);
        thrust::copy(pts_y.begin(), pts_y.begin() + cnp, host_pts_y.begin() + pnp);
        thrust::copy(pts_z.begin(), pts_z.begin() + cnp, host_pts_z.begin() + pnp);
        thrust::copy(nrs_x.begin(), nrs_x.begin() + cnp, host_nrs_x.begin() + pnp);
        thrust::copy(nrs_y.begin(), nrs_y.begin() + cnp, host_nrs_y.begin() + pnp);
        thrust::copy(nrs_z.begin(), nrs_z.begin() + cnp, host_nrs_z.begin() + pnp);

        processed += p;

        std::cout << "num processed: " << processed << std::endl;
    }

    run_time_points_and_normals.stop();

    // Free up memory
    image_data.resize(0);

    pts_x.resize(0);   nrs_x.resize(0);
    pts_y.resize(0);   nrs_y.resize(0);
    pts_z.resize(0);   nrs_z.resize(0);

    ax.resize(0);    bx.resize(0);    cx.resize(0);
    ay.resize(0);    by.resize(0);    cy.resize(0);
    az.resize(0);    bz.resize(0);    cz.resize(0);
    axn.resize(0);   bxn.resize(0);   cxn.resize(0);
    ayn.resize(0);   byn.resize(0);   cyn.resize(0);
    azn.resize(0);   bzn.resize(0);   czn.resize(0);

    ///////////////////////////////////////////////////////////////////////////
    // Set and allocate triangles

    util::Timer run_time_triangles;

    vector<int> trs0(num_triangles);
    vector<int> trs1(num_triangles);
    vector<int> trs2(num_triangles);

    for_each(
        policy,
        make_counting_iterator(0),
        make_counting_iterator((nx-1)*(ny-1)*(nz-1)),
        set_triangles(
            trs0.data(), trs1.data(), trs2.data(),
            nx, ny, nz,
            a0.data(), b0.data(), c0.data(),
            cube_ids.data(),
            tri_scan.data()));

    host_vector<int> host_trs_0(trs0);
    host_vector<int> host_trs_1(trs1);
    host_vector<int> host_trs_2(trs2);

    run_time_triangles.stop();
    run_time.stop();

    doc.add("Number of vertices in mesh", num_points);
    doc.add("Number of triangles in mesh", num_triangles);

    doc.add("Allocate Memory", "");
    doc.get("Allocate Memory")->add(
        "CPU Time (clicks)", run_time_allocate_memory.getTotalTicks());
    doc.get("Allocate Memory")->add(
        "CPU Time (seconds)", run_time_allocate_memory.getCPUtime());
    doc.get("Allocate Memory")->add(
        "Wall Time (seconds)", run_time_allocate_memory.getWallTime());

    doc.add("Set Points and Normals", "");
    doc.get("Set Points and Normals")->add(
        "CPU Time (clicks)", run_time_points_and_normals.getTotalTicks());
    doc.get("Set Points and Normals")->add(
        "CPU Time (seconds)", run_time_points_and_normals.getCPUtime());
    doc.get("Set Points and Normals")->add(
        "Wall Time (seconds)", run_time_points_and_normals.getWallTime());

    doc.add("Set Triangles", "");
    doc.get("Set Triangles")->add(
        "CPU Time (clicks)", run_time_triangles.getTotalTicks());
    doc.get("Set Triangles")->add(
        "CPU Time (seconds)", run_time_triangles.getCPUtime());
    doc.get("Set Triangles")->add(
        "Wall Time (seconds)", run_time_triangles.getWallTime());

    doc.add("Total Program CPU Time (clicks)", run_time.getTotalTicks());
    doc.add("Total Program CPU Time (seconds)", run_time.getCPUtime());
    doc.add("Total Program WALL Time (seconds)", run_time.getWallTime());

    std::cout << doc.generateYAML();

    util::saveTriangleMesh(
        outFile,
        host_pts_x, host_pts_y, host_pts_z,
        host_nrs_x, host_nrs_y, host_nrs_z,
        host_trs_0, host_trs_1, host_trs_2);

/*
 ******* Metrics to test correctness*******
    scalar_t sum_pts = 0.0;
    scalar_t sum_nrs = 0.0;
    for(int w = 0; w != num_points; ++w)
    {
        sum_pts += host_pts_x[w] + host_pts_y[w] + host_pts_z[w];
        sum_nrs += host_nrs_x[w] + host_nrs_y[w] + host_nrs_z[w];

        while(sum_pts > 500000000)
            sum_pts -=  500000000;
        while(sum_pts <-500000000)
            sum_pts +=  500000000;
        while(sum_nrs > 500000000)
            sum_nrs -=  500000000;
        while(sum_nrs <-500000000)
            sum_nrs +=  500000000;
    }

    scalar_t sum_tri_pts = 0.0;
    for(int v = 0; v != num_triangles; ++v)
    {
        int p1 = host_trs0[v];
        int p2 = host_trs1[v];
        int p3 = host_trs2[v];

        sum_tri_pts += host_pts_x[p1] + host_pts_y[p1] + host_pts_z[p1] +
                       host_pts_x[p2] + host_pts_y[p2] + host_pts_z[p2] +
                       host_pts_x[p3] + host_pts_y[p3] + host_pts_z[p3];


        while(sum_tri_pts > 500000000)
            sum_tri_pts -=  500000000;
        while(sum_tri_pts <-500000000)
            sum_tri_pts +=  500000000;
    }

    std::cout << "num points " << num_points << ", " << host_pts_x.size() << std::endl;
    std::cout << "num triangles " << num_triangles << ", " << host_trs0.size() << std::endl;
    std::cout << "sum points " << sum_pts << std::endl;
    std::cout << "sum normals " << sum_nrs << std::endl;
    std::cout << "sum triangle points " << sum_tri_pts << std::endl;
*/
}


