#include <iostream>
#include <string.h>
#include <cstdlib>

#include "FlyingEdgesAlgorithm.h"

#include "../util/LoadImage.h"
#include "../util/SaveTriangleMesh.h"

#include "../util/Timer.h"
#include "../mantevoCommon/YAML_Doc.hpp"

__global__ void doubler(scalar_t* a)
{
    int i = blockIdx.x;

    a[i] = a[i]*2;
}

int main(int argc, char* argv[])
{
/*
    std::vector<scalar_t> vec(10);
    scalar_t c = 9.3;
    for(scalar_t& s: vec)
    {
        s = c;
        c += 1.15;
    }

    scalar_t* deviceVec;
    cudaMalloc(&deviceVec, 10*sizeof(scalar_t));

    cudaMemcpy(deviceVec, vec.data(), 10*sizeof(scalar_t), cudaMemcpyHostToDevice);

    doubler<<<10,1>>>(deviceVec);

    cudaMemcpy(vec.data(), deviceVec, 10*sizeof(scalar_t), cudaMemcpyDeviceToHost);

    cudaFree(deviceVec);

    for(scalar_t const& s: vec)
        std::cout << s << std::endl;

    return 0;
*/

    float isoval;
    bool isovalSet = false;
    char* vtkFile = NULL;
    char* outFile = NULL;
    std::string yamlDirectory = "";
    std::string yamlFileName  = "";

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

    YAML_Doc doc("Flying Edges", "0.1", yamlDirectory, yamlFileName);

    doc.add("Flying Edges Algorithm", "cuda");
    doc.add("Volume image data file path", vtkFile);
    doc.add("Polygonal mesh output file", outFile);
    doc.add("Isoval", isoval);

    util::Image3D image = util::loadImage(vtkFile);

    doc.add("File x-dimension", image.xdimension());
    doc.add("File y-dimension", image.ydimension());
    doc.add("File z-dimension", image.zdimension());

    util::Timer runTime;

    util::Timer toDeviceTime;
    FlyingEdgesAlgorithm algo(image, isoval);
    toDeviceTime.stop();

    util::Timer runTimePass1;
    algo.pass1();
    runTimePass1.stop();

    util::Timer runTimePass2;
    algo.pass2();
    runTimePass2.stop();

    util::Timer runTimePass3;
    algo.pass3();
    runTimePass3.stop();

    util::Timer runTimePass4;
//    algo.pass4();
    runTimePass4.stop();

    util::Timer fromDeviceTime;
    util::TriangleMesh mesh = algo.moveOutput();
    fromDeviceTime.stop();

    runTime.stop();

    doc.add("Number of vertices in mesh", mesh.numberOfVertices());
    doc.add("Number of triangles in mesh", mesh.numberOfTriangles());

    doc.add("To Device", "");
    doc.get("To Device")->add("CPU Time (clicks)", toDeviceTime.getTotalTicks());
    doc.get("To Device")->add("CPU Time (seconds)", toDeviceTime.getCPUtime());
    doc.get("To Device")->add("Wall Time (seconds)", toDeviceTime.getWallTime());

    doc.add("Pass 1", "");
    doc.get("Pass 1")->add("CPU Time (clicks)", runTimePass1.getTotalTicks());
    doc.get("Pass 1")->add("CPU Time (seconds)", runTimePass1.getCPUtime());
    doc.get("Pass 1")->add("Wall Time (seconds)", runTimePass1.getWallTime());

    doc.add("Pass 2", "");
    doc.get("Pass 2")->add("CPU Time (clicks)", runTimePass2.getTotalTicks());
    doc.get("Pass 2")->add("CPU Time (seconds)", runTimePass2.getCPUtime());
    doc.get("Pass 2")->add("Wall Time (seconds)", runTimePass2.getWallTime());

    doc.add("Pass 3", "");
    doc.get("Pass 3")->add("CPU Time (clicks)", runTimePass3.getTotalTicks());
    doc.get("Pass 3")->add("CPU Time (seconds)", runTimePass3.getCPUtime());
    doc.get("Pass 3")->add("Wall Time (seconds)", runTimePass3.getWallTime());

    doc.add("Pass 4", "");
    doc.get("Pass 4")->add("CPU Time (clicks)", runTimePass4.getTotalTicks());
    doc.get("Pass 4")->add("CPU Time (seconds)", runTimePass4.getCPUtime());
    doc.get("Pass 4")->add("Wall Time (seconds)", runTimePass4.getWallTime());

    doc.add("From Device", "");
    doc.get("From Device")->add("CPU Time (clicks)", fromDeviceTime.getTotalTicks());
    doc.get("From Device")->add("CPU Time (seconds)", fromDeviceTime.getCPUtime());
    doc.get("From Device")->add("Wall Time (seconds)", fromDeviceTime.getWallTime());

    doc.add("Total Program CPU Time (clicks)", runTime.getTotalTicks());
    doc.add("Total Program CPU Time (seconds)", runTime.getCPUtime());
    doc.add("Total Program WALL Time (seconds)", runTime.getWallTime());

    std::cout << doc.generateYAML();

//    util::saveTriangleMesh(mesh, outFile);
}
