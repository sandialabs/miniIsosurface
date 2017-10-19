/*
 * vtkReference/main.cpp
 *
 *  Created on: Jun 28, 2017
 *      Author: dbourge
 *
 * miniIsosurface is distributed under the OSI-approved BSD 3-clause License.
 * See LICENSE.txt for details.
 *
 * Copyright (c) 2017
 * National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under
 * the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
 * certain rights in this software.
 */
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkMarchingCubes.h>
#include <vtkFlyingEdges3D.h>
#include <vtkPolyDataWriter.h>

#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPoints.h>
#include <vtkImageData.h>
#include <vtkDataSet.h>
#include <vtkDataObject.h>
#include <vtkExtractVOI.h>

#include <iostream>

#include "../marchingCubes/util/Timer.h"
#include "../marchingCubes/mantevoCommon/YAML_Doc.hpp"

int main(int argc, char *argv[])
{
    float isoval;
    bool isovalSet = false;
    char* vtkFile = NULL;
    char* outFile = NULL;
    std::string yamlDirectory = "";
    std::string yamlFileName  = "";
    bool useMC = true;

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
        if( (strcmp(argv[i], "-a") == 0) || (strcmp(argv[i], "-algorithm") == 0))
        {
            useMC = ("mc" == argv[++i]);
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
                "Serial Marching Cubes Options:"  << std::endl <<
                "  -input_file (-i)"              << std::endl <<
                "  -input_dat"                    << std::endl <<
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
    YAML_Doc doc("VTK Reference", "0.2", yamlDirectory, yamlFileName);

    // Add information related to this run to doc.
    std::string algo = useMC ? "Marching Cubes" : "Flying Edges";
    doc.add("VTK Reference", algo);
    doc.add("Volume image data file path", vtkFile);
    doc.add("Polygonal mesh output file", outFile);
    doc.add("Isoval", isoval);

    // Load the image file
    vtkSmartPointer<vtkStructuredPointsReader> reader =
        vtkSmartPointer<vtkStructuredPointsReader>::New();
    reader->SetFileName(vtkFile);
    reader->Update();

    vtkSmartPointer<vtkImageData> volume =
        vtkSmartPointer<vtkImageData>::New();
    volume->DeepCopy(reader->GetOutput());

    int dims[3];
    volume->GetDimensions(dims);

    doc.add("File x-dimension", dims[0]);
    doc.add("File y-dimension", dims[1]);
    doc.add("File z-dimension", dims[2]);

    // Time the output. Timer's constructor starts timing.
    util::Timer runTime;

    // Run the algorithm
    vtkSmartPointer<vtkPolyDataWriter> writer =
        vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName(outFile);
    vtkPolyData* polygonalMesh;

    if(useMC)
    {
        vtkSmartPointer<vtkMarchingCubes> surface =
            vtkSmartPointer<vtkMarchingCubes>::New();

#if VTK_MAJOR_VERSION <= 5
        surface->SetInput(volume);
#else
        surface->SetInputData(volume);
#endif
        surface->ComputeNormalsOn();
        surface->SetValue(0, isoval);
        surface->Update();

        writer->SetInputConnection(surface->GetOutputPort());
        writer->Update();
        polygonalMesh = surface->GetOutput(0);
    }
    else
    {
        vtkSmartPointer<vtkFlyingEdges3D> surface =
            vtkSmartPointer<vtkFlyingEdges3D>::New();

#if VTK_MAJOR_VERSION <= 5
        surface->SetInput(volume);
#else
        surface->SetInputData(volume);
#endif
        surface->ComputeNormalsOn();
        surface->SetValue(0, isoval);
        surface->Update();

        writer->SetInputConnection(surface->GetOutputPort());
        writer->Update();
        polygonalMesh = surface->GetOutput(0);
    }

    // End timing
    runTime.stop();

    // Report mesh information
    doc.add("Number of triangles in mesh", polygonalMesh->GetNumberOfPolys());

    // Report timing information
    doc.add("Total Program CPU Time (clicks)", runTime.getTotalTicks());
    doc.add("Total Program CPU Time (seconds)", runTime.getCPUtime());
    doc.add("Total Program WALL Time (seconds)", runTime.getWallTime());

    // Generate the YAML file. The file will be both saved and printed to console.
    std::cout << doc.generateYAML();

    // Save the polygonal mesh to the output file
    writer->Write();
}
