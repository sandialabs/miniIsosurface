/*
 * serial/main.cpp
 *
 *  Created on: Feb 20, 2017
 *      Author: dbourge
 */

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include "open-simplex-noise.h"

#include "../../flyingEdges/util/ConvertBuffer.h"

#include <iostream>
#include <fstream>
#include <utility>
#include <array>
#include <vector>
#include <string.h>

std::pair<std::vector<float>, std::array<float, 3> > dataGen(
    int ni, int nj, int nk,
    float mask_thres, float bot_mask_thres,
    double noisespacefreq, double noisetimefreq);

void writeFile(
    const char* outFile,
    std::vector<float> data,
    std::array<float, 3> spacing,
    int ni, int nj, int nk);

int main(int argc, char* argv[])
{
    char* outFile = (char*)"dataGen.vtk";

    int ni = 0;                   /* Global grid size */
    int nj = 0;
    int nk = 0;

    float mask_thres=0.0;     /* upper mask threshold  (range -1 to 1) */
    float bot_mask_thres=0.5; /* bottom mask threshold (range 0.0 to (mask_thres+1)/2 ) */

    double noisespacefreq = 10.0; /* Spatial frequency of noise */
    double noisetimefreq = 0.25;  /* Temporal frequency of noise */

    int tstart = 0;

    // Read command line arguments
    for(int i=0; i<argc; i++)
    {
        if( (strcmp(argv[i], "-o") == 0) || (strcmp(argv[i], "-outfile") == 0))
        {
            outFile = argv[++i];
        }
        else if( (strcmp(argv[i], "-size") == 0))
        {
            ni = atoi(argv[++i]);
            nj = atoi(argv[++i]);
            nk = atoi(argv[++i]);
        }
        else if( (strcmp(argv[i], "-m") == 0) || (strcmp(argv[i], "-maskthreshold") == 0))
        {
            mask_thres = strtof(argv[++i], NULL);
            bot_mask_thres = (mask_thres+1.0)/2;
        }
        else if( (strcmp(argv[i], "-s") == 0) || (strcmp(argv[i], "-noisespacefre") == 0))
        {
            noisespacefreq = strtod(argv[++i], NULL);
        }
        else if( (strcmp(argv[i], "-t") == 0) || (strcmp(argv[i], "-noisetimefreq") == 0))
        {
            noisetimefreq = strtod(argv[++i], NULL);
        }
        else if( (strcmp(argv[i], "-a") == 0) || (strcmp(argv[i], "-tstart") == 0))
        {
            tstart = atoi(argv[++i]);
        }
        else if( (strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "-help") == 0))
        {
            std::cout <<
                "Usage: ./dataGen -size NI NJ NK [options]"     << std::endl <<
                "DataGen Options: "                             << std::endl <<
                "  -outfile (-o)"                               << std::endl <<
                "   write to this file, default dataGen.vtk"    << std::endl <<
                "  -maskthreshold (-m)"                         << std::endl <<
                "   valid values are floats between -1.0 and "
                    "1.0, default 0.0"                          << std::endl <<
                "  -noisespacefreq (-s)"                        << std::endl <<
                "   spatial frequency of noise function "
                    "default 10.0"                              << std::endl <<
                "  -noisetimefreq (-t)"                         << std::endl <<
                "   spatial frequency of noise function "       // does anything?
                    "default 0.25"                              << std::endl <<
                "  -tstart (-a)"                                << std::endl <<
                "   starting time step, valid values are >= 0 "
                    "default 0"                                 << std::endl;
            return 0;
        }
    }

    if(ni == 0 || nj == 0 || nk == 0)
    {
        std::cout <<
            "Error. Usage: ./dataGen -size NI NJ NK [options]" << std::endl <<
            "Try -help" << std::endl;
        return 0;
    }

    std::cout << "Generating data" << std::endl;

    std::pair<std::vector<float>, std::array<float, 3> > data = dataGen(
            ni, nj, nk,
            mask_thres, bot_mask_thres,
            noisespacefreq, noisetimefreq);

    std::cout << "Generated data" << std::endl;

//    int count = 0;
//    for(auto val: data.first)
//        if(val != -999)
//            count += 1;
//    std::cout << count << " / " << ni*nj*nk << std::endl;

    writeFile(outFile, std::move(data.first), data.second, ni, nj, nk);
}

// This function is lifted from miniIO/struct.c and changed for purposes of this
// miniapp.
std::pair<std::vector<float>, std::array<float, 3> >
dataGen(
    int ni, int nj, int nk,
    float mask_thres, float bot_mask_thres,
    double noisespacefreq, double noisetimefreq)
{
    std::vector<float> ret;
    int a, i, j, k, t;            /* loop indices */
    int tt;                       /* Actual time step from tstart */
    float x, y, z;
    int tstart = 0;
    int nt = 1;                  /* Number of time steps */
    int inp = 1;      /* Number of tasks in i */
    int jnp = 1;      /* Number of tasks in j */
    int knp = 1;      /* Number of tasks in k */
    int numtask = 1;               /* task per processor */
    int point_id, tmp_id;          /* grid point id */
    int x_index, y_index, z_index; /* point index along each axis */
    int xy_dims,  x_dims;
    float deltax, deltay, deltaz;
    int numPoints;
    float *data;
    float *height;
    int height_index;
    int hindex;
    int maskTindex;
    int *ola_mask;
    int *ol_mask;
    int mask_thres_index;
    struct osn_context *simpnoise;    /* Open simplex noise context */

    const int num_varnames=4;
    char *varnames[num_varnames];

    /* MPI vars */
    int cprocs[3], cpers[3], crnk[3];  /* MPI Cartesian info */
    int rank = 0, nprocs = 1;
    int cni, cnj, cnk;   /* Points in this task */
    int is, js, ks;     /* Global index starting points */
    float xs, ys, zs;    /* Global coordinate starting points */

    numPoints = ni*nj*nk;

    ret.resize(numPoints);
    data = ret.data();

     /* Set up Cartesian communicator */
    cprocs[0] = inp;  cprocs[1] = jnp;  cprocs[2] = knp;
    cpers[0] = 0;  cpers[1] = 0;  cpers[2] = 0;    /* No periodicity */

    deltax = 1.f/(ni-1);
    deltay = 1.f/(nj-1);
    deltaz = 1.f/(nk-1);
    cni = ni / inp;
    cnj = nj / jnp;
    cnk = nk / knp;
    is = crnk[0] * cni;
    js = crnk[1] * cnj;
    ks = crnk[2] * cnk;
    xs = is * deltax;
    ys = js * deltay;
    zs = ks * deltaz;

    xy_dims = ni * nj;
    x_dims = ni;

    /* adjust mask threshold  to compensate by bottom threshold */
    mask_thres = mask_thres - bot_mask_thres;
    mask_thres_index = (int) ( ((mask_thres+1)/2) * (nk-1));
    maskTindex = nk-1;

    /* Set up osn */
    open_simplex_noise(1234, &simpnoise);   /* Fixed seed, for now */

    /* Allocate arrays */
    height = (float *) malloc((size_t)cni*cnj*cnk*sizeof(float));
    ola_mask = (int *) malloc((size_t)cni*cnj*cnk*sizeof(int));
    ol_mask = (int *) malloc((size_t)cni*cnj*cnk*sizeof(int));

    //varnames[0] = "data";
    //varnames[1] = "height";
    //varnames[2] = "ola_mask";
    //varnames[3] = "ol_mask";

    /* generate masked grid */
    /* Spatial loops */
    size_t ii;     /* data index */

    z = zs;
    for(k = 0, ii = 0; k < cnk; k++) {
        y = ys;
        for(j = 0; j < cnj; j++) {
            x = xs;
            for(i = 0; i < cni; i++, ii++) {
                x_index = (int) (x/deltax);
                y_index = (int) (y/deltay);
                z_index = (int) (z/deltaz);

                /* calculate point index */
                point_id = (z_index * xy_dims) + (y_index * x_dims) + x_index;

                /* Get height and subtract bottom threshold */
                height[ii] =  (float)open_simplex_noise2(simpnoise, x*noisespacefreq, y*noisespacefreq)  - bot_mask_thres;

                /* height_index = (int) height[ii]/deltaz; */
                height_index = (int) (((height[ii]+1)/2) * (nk-1));

                /* Calculate ola_mask values */
                if (z_index > mask_thres_index  && z_index > height_index) {
                    ola_mask[ii] = 2;  /* Atmosphere */
                }
                else if (z_index < mask_thres_index  && z_index > height_index) {
                    ola_mask[ii] = 0;  /* ocean */
                }
                else if (z_index <= height_index) {
                    if (height[ii] >= mask_thres  || z_index < height_index) {
                        ola_mask[ii] = 1;  /* land */
                    } else {
                        ola_mask[ii] = 0;  /* ocean */
                    }
                }
                else if (z_index == mask_thres_index  && height[ii] <= mask_thres) {
                        ola_mask[ii] = 0;  /* ocean */
                    }
                else {
                      printf("WARNING: ola_mask condition not considered for Point_index: (%d,%d,%d)\n"
                     "Point_id: %d  Height: %f HeightID: %d  mask_thres_index=%d\n",
                     x_index, y_index, z_index, point_id+1, height[ii], height_index, mask_thres_index);
                }


                hindex = (int) ( (((height[ii]+1)/2) * (nk-1)) / (mask_thres_index+1)) * (nk-1);

                hindex = (int) ((height[ii]+1) / (mask_thres+1) * (nk-1));

                if (hindex > maskTindex) {
                    hindex = maskTindex;
                }

                /* Calculate ol_mask values */
                if (z_index < maskTindex  && z_index > hindex) {
                    ol_mask[ii] = 0;  /* ocean */
                }
                else if (z_index <= hindex) {
                    if (height[ii] >= mask_thres  || z_index < hindex) {
                        ol_mask[ii] = 1;  /* land */
                     }
                else {
                    ol_mask[ii] = 0;  /* ocean */
                  }
                }
                else if (z_index == maskTindex  && height[ii] <= mask_thres) {
                ol_mask[ii] = 0;  /* ocean */
                }
                else {
                     printf("WARNING: ol_mask condition not considered for Point_index: (%d,%d,%d)\n"
                     "Point_id: %d  Height: %f HeightID: %d maskTindex=%d\n",
                     x_index, y_index, z_index, point_id+1, height[ii], hindex, maskTindex);
                }

                x += deltax;
              }
              y += deltay;
            }
            z += deltaz;
          }

          /* generate ocean land data */
      for(t = 0, tt = tstart; t < nt; t++, tt++) {
        /* Spatial loops */

        z = zs;
        for(k = 0, ii = 0; k < cnk; k++) {
          y = ys;
          for(j = 0; j < cnj; j++) {
        x = xs;
            for(i = 0; i < cni; i++, ii++) {

          if ( ol_mask[ii] == 0) {
            /* if  ( ola_mask[ii] == 0) { */
            data[ii] = (float)open_simplex_noise4(simpnoise, x*noisespacefreq, y*noisespacefreq, z*noisespacefreq, tt*noisetimefreq);
          }
          else {
            data[ii] = -999; // FILLVALUE
          }

          x += deltax;
        }
        y += deltay;
          }
          z += deltaz;
        }





    open_simplex_noise_free(simpnoise);
    free(height);
    free(ola_mask);
    free(ol_mask);
}

    return std::make_pair(ret, std::array<float, 3>({deltax, deltay, deltaz}));
}

void writeFile(
    const char* outFile,
    std::vector<float> data,
    std::array<float, 3> spacing,
    int nx, int ny, int nz)
{
    std::cout << "Writing file..." << std::endl;

    std::ofstream outStream(outFile);

    outStream << "# vtk DataFile Version 3.0" << std::endl;
    outStream << "vtk output" << std::endl;
    outStream << "BINARY" << std::endl;
    outStream << "DATASET STRUCTURED_POINTS" << std::endl;
    outStream << "DIMENSIONS " << nx << " " << ny << " " << nz << std::endl;
    outStream << "SPACING " << spacing[0] << " " << spacing[1] << " "
           << spacing[2] << std::endl;
    outStream << "ORIGIN " << 0 << " " << 0 << " "
           << 0 << std::endl;
    outStream << "POINT_DATA " << nx*ny*nz << std::endl;
    outStream << "SCALARS ImageFile float" << std::endl;
    outStream << "LOOKUP_TABLE default" << std::endl;

    std::vector<char> wbuff;
    size_t step = 5000000;
    for(size_t i = 0; i < data.size(); i += step)
    {
        std::cout << i << " / " << data.size() << std::endl;

        size_t bufsize = std::min(data.size() - i, step);
        bufsize = bufsize * 3 * sizeof(float);
        wbuff.resize(bufsize);
        float* bufPointer = reinterpret_cast<float*>(wbuff.data());

        for(int idx = i; idx != std::min(data.size(), i + step); ++idx)
        {
            *bufPointer = data[idx];
            util::flipEndianness(*bufPointer++);
        }

        outStream.write(wbuff.data() + i, wbuff.size());
    }
    outStream << std::endl;
    outStream.close();

//    std::ofstream stream(outFile);
//
//    stream << "# vtk DataFile Version 3.0" << std::endl;
//    stream << "vtk output" << std::endl;
//    stream << "BINARY" << std::endl;
//    stream << "DATASET STRUCTURED_POINTS" << std::endl;
//    stream << "DIMENSIONS " << ni << " " << nj << " " << nk << std::endl;
//    stream << "SPACING " << spacing[0] << " " << spacing[1] << " "
//           << spacing[2] << std::endl;
//    stream << "ORIGIN " << 0 << " " << 0 << " " << 0 << std::endl;
//    stream << "POINT_DATA " << ni*nj*nk << std::endl;
//    stream << "SCALARS ImageFile float" << std::endl;
//    stream << "LOOKUP_TABLE default" << std::endl;
//
//    std::vector<char> wbuff;
//    std::size_t bufsize = ni*nj*nk*3*sizeof(float);
//    wbuff.resize(bufsize);
//
//    float* bufPointer = reinterpret_cast<float*>(wbuff.data());
//
//    for(float const& val: data)
//    {
//        *bufPointer = val;
//        util::flipEndianness(*bufPointer++);
//    }
//    stream.write(wbuff.data(), wbuff.size());
//    stream << std::endl;
//    stream.close();
}































