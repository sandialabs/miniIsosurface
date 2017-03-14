#include <iostream>
#include <string.h>
#include <cstdlib>

#include "FlyingEdgesAlgorithm.h"

#include "../util/LoadImage.h"
#include "../util/SaveTriangleMesh.h"

#include "../util/Timer.h"
#include "../mantevoCommon/YAML_Doc.hpp"

int main(int argc, char* argv[])
{
    scalar_t isoval;
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

    // Create a yamlDoc. If yamlDirectory and yamlFileName weren't assigned,
    // YAML_Doc will create a file at in the current directory with a
    // timestamp on it.
    YAML_Doc doc("Flying Edges", "0.1", yamlDirectory, yamlFileName);

    // Add information related to this run to doc.
    doc.add("Flying Edges Algorithm", "serial");
    doc.add("Volume image data file path", vtkFile);
    doc.add("Polygonal mesh output file", outFile);
    doc.add("Isoval", isoval);

    // Load the image file
    util::Image3D image = util::loadImage(vtkFile);

    doc.add("File x-dimension", image.xdimension());
    doc.add("File y-dimension", image.ydimension());
    doc.add("File z-dimension", image.zdimension());

    // Time the output. util::Timer's constructor starts timing.
    util::Timer runTime;

    // The inputs of the algorithm are the 3D image file and the isoval to
    // estimate the isosurface at.
    FlyingEdgesAlgorithm algo(image, isoval);
    // The flying edges algorithm makes 4 passes through the image file.
    // Each pass is timed.

    // Pass 1 of the algorithm labels each edge parallel to the x-axis as cut
    // or not. In the process, each gridEdge is assigned an xl and xr.
    // All edges before xl are uncut and all edges after xr are uncut.
    // Subsequent passes of the algorithm don't look outside of [xl, xr).
    // A gridEdge E_jk can be thought of as the row of edges parallel to the
    // x-axis for some fixed j and k.
    util::Timer runTimePass1;
    algo.pass1();
    runTimePass1.stop();

    // Pass 2 of the algorithm determines the marching cubes case ID of each
    // cube. This is determined fully from information obtained in pass1, so
    // there is no need to access the input image. Each cube starts at (i,j,k)
    // and extends to (i+1, j+1, k+1).
    // In addition to determining case ID of each cell, pass 2 counts the
    // number of cuts on incident to each gridEdge.
    util::Timer runTimePass2;
    algo.pass2();
    runTimePass2.stop();

    // Pass 3 of the algorithm uses information from pass 2 to determine how
    // many triangles and points there are. It also sets up starting indices
    // on each gridEdge. Once these sizes are determined, memory is allocated
    // for storing triangles, points and normals.
    util::Timer runTimePass3;
    algo.pass3();
    runTimePass3.stop();

    // Pass 4 of the algorithm calculates calculates and fills out points,
    // normals and the triangles.
    util::Timer runTimePass4;
    algo.pass4();
    runTimePass4.stop();

    // This function receives the output. The data is not copied or deep copied
    // but instead is moved or shallow copied. Once moveOutput is called, the
    // algo structure no longer maintains responsibility of any data.
    util::TriangleMesh mesh = algo.moveOutput();

    // End overall timing
    runTime.stop();

    // Report mesh information
    doc.add("Number of vertices in mesh", mesh.numberOfVertices());
    doc.add("Number of triangles in mesh", mesh.numberOfTriangles());

    // Report timing information
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

    doc.add("Total Program CPU Time (clicks)", runTime.getTotalTicks());
    doc.add("Total Program CPU Time (seconds)", runTime.getCPUtime());
    doc.add("Total Program WALL Time (seconds)", runTime.getWallTime());

    // Generate the YAML file. The file will be both saved and printed to console.
    std::cout << doc.generateYAML();

    // Save the polygonal mesh to the output file.
//    util::saveTriangleMesh(mesh, outFile);
}
