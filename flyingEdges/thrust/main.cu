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
// Next: Split the image up and contain ghost sell information.
//
// The problem is that for large input, the entire image will not fit
// into memory.
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

    // tri_scan
    vector<int> tri_scan(n);
    vector<uchar> cube_ids(n);

    vector<int> a0(n);
    vector<int> b0(n);
    vector<int> c0(n);

    int num_triangles = 0;
    int num_points = 0;

    int max_cur_points = process_size / 100; // used for keeping track of amount
                                             // reserved in pts, nrs vectors
    vector<scalar_t> pts_x(max_cur_points);  // process_size is too big of a guess because
    vector<scalar_t> pts_y(max_cur_points);  // that means every single (x,y,z) is cut.
    vector<scalar_t> pts_z(max_cur_points);  // going with process_size/100 for now.

    vector<scalar_t> nrs_x(max_cur_points);
    vector<scalar_t> nrs_y(max_cur_points);
    vector<scalar_t> nrs_z(max_cur_points);

    int guess_num_pts = n / 100;            // Again, just a guess for amount to reserve.
    host_vector<scalar_t> host_pts_x;     host_pts_x.reserve(guess_num_pts);
    host_vector<scalar_t> host_pts_y;     host_pts_y.reserve(guess_num_pts);
    host_vector<scalar_t> host_pts_z;     host_pts_z.reserve(guess_num_pts);
    host_vector<scalar_t> host_nrs_x;     host_nrs_x.reserve(guess_num_pts);
    host_vector<scalar_t> host_nrs_y;     host_nrs_y.reserve(guess_num_pts);
    host_vector<scalar_t> host_nrs_z;     host_nrs_z.reserve(guess_num_pts);

    run_time_allocate_memory.stop();

    // ----- The algorihtm in a nutshel -----
    // While not all the image has been processed:
    //   First:  fill out pre scan values. That is a0, b0, c0, tri_scan, case_id
    //   Second: scan a0, b0, c0, tri_scan.
    //   Third:  Calculate points and normals. Put to output data at a0, b0 and c0 vals.
    // Free up image data
    // Calculate triangles from a0, b0, c0, tri_scan and case_id
    // --------------------------------------

    util::Timer run_time_points_and_normals;

    std::cout << "num to process: " << n << std::endl;

    while(processed != n)
    {
        int p = std::min(process_size, n - processed);

        // calculate a0, b0, c0, tri_scan

        auto scan_iterator = make_zip_iterator(
            make_tuple(
                a0.begin() + processed,
                b0.begin() + processed,
                c0.begin() + processed,
                tri_scan.begin() + processed,
                cube_ids.begin() + processed));

        // pass1 calculate whether or not a0, b0, c0 is cut, the number of triangles
        // at the (x,y,z) cube and the cube_id.

        transform(
            policy,
            make_counting_iterator(processed),
            make_counting_iterator(processed + p),
            scan_iterator,
            fill_out_pre_scan_values(
                nx, ny, nz,
                isoval,
                image_data.data()));

        ///////////////////////////////////////////////////////////////////////
        // pass 2: scan step + allocate for points and normal values
        // tmp sums
        int tmp_a = a0[processed + p - 1];
        int tmp_b = b0[processed + p - 1];
        int tmp_c = c0[processed + p - 1];
        int tmp_t = tri_scan[processed + p - 1];

        exclusive_scan(
            policy,
            a0.begin() + processed,
            a0.begin() + processed + p,
            a0.begin() + processed,
            num_points);                    // increase starting value

        exclusive_scan(
            policy,
            b0.begin() + processed,
            b0.begin() + processed + p,
            b0.begin() + processed,
            tmp_a + a0[processed + p - 1]); // increase starting value

        exclusive_scan(
            policy,
            c0.begin() + processed,
            c0.begin() + processed + p,
            c0.begin() + processed,
            tmp_b + b0[processed + p - 1]); // increase starting value

        exclusive_scan(
            policy,
            tri_scan.begin() + processed,
            tri_scan.begin() + processed + p,
            tri_scan.begin() + processed,
            num_triangles);

        int prev_num_points = num_points;

        num_triangles = tmp_t + tri_scan[processed + p - 1];
        num_points = tmp_c + c0[processed + p - 1];

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

        ///////////////////////////////////////////////////////////////////////
        // pass 3: calculate points and normal values.
        // Triangles will be done at the end.

        auto pts_nors_info_beg =
            make_zip_iterator(
                make_tuple(
                    make_counting_iterator(processed),
                    scan_iterator));
        auto pts_nors_info_end = pts_nors_info_beg + p;

        for_each(
            policy,
            pts_nors_info_beg,
            pts_nors_info_end,
            calculate_points_and_normals(
                nx, ny, nz,
                spacing_x, spacing_y, spacing_z,
                zeropos_x, zeropos_y, zeropos_z,
                isoval,
                image_data.data(),
                prev_num_points,
                pts_x.data(), pts_y.data(), pts_z.data(),
                nrs_x.data(), nrs_y.data(), nrs_z.data()));

        // TODO make sure host_pts_x has a good initial guess to the size.
        // Its on the host, so can make big lar
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

        ///////////////////////////////////////////////////////////////////////
        processed += p;
        std::cout << "num processed: " << processed << std::endl;
    }

    run_time_points_and_normals.stop();

    // Free up memory
    image_data.resize(0);

    pts_x.resize(0);   nrs_x.resize(0);
    pts_y.resize(0);   nrs_y.resize(0);
    pts_z.resize(0);   nrs_z.resize(0);

    ///////////////////////////////////////////////////////////////////////////
    // Set and allocate triangles
    util::Timer run_time_set_triangles;

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

    run_time_set_triangles.stop();
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
        "CPU Time (clicks)", run_time_set_triangles.getTotalTicks());
    doc.get("Set Triangles")->add(
        "CPU Time (seconds)", run_time_set_triangles.getCPUtime());
    doc.get("Set Triangles")->add(
        "Wall Time (seconds)", run_time_set_triangles.getWallTime());

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


