///////////////////////////////////////////////////////////////////////////////
// *Time: 5e-5 seconds
///////////////////////////////////////////////////////////////////////////////
//

/*
__global__
void pass1gpu(
    scalar_t* pointValues,  // input
    int nx, int ny, int nz, // input
    scalar_t isoval,        // input
    uchar* edgeCases)       // output
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if(j >= ny || k >= nz)
        return;

    edgeCases[k*ny*(nx-1) + j*(nx-1) + i] =
        calcCaseEdge(
            pointValues[k*ny*nx + j*nx + i + 0] >= isoval,
            pointValues[k*ny*nx + j*nx + i + 1] >= isoval);
}
*/

/*
__global__
void pass1gpu(
    scalar_t* pointValues,  // input
    scalar_t isoval,        // input
    uchar* edgeCases)       // output
{
    int i = blockIdx.x;
    int j = blockIdx.y;
    int k = blockIdx.z;

    int ny = blockDim.y;
    int nx = blockDim.x + 1;

    edgeCases[k*ny*(nx-1) + j*(nx-1) + i] =
        calcCaseEdge(
            pointValues[k*ny*nx + j*nx + i + 0] >= isoval,
            pointValues[k*ny*nx + j*nx + i + 1] >= isoval);
}

void FlyingEdgesAlgorithm::pass1()
{
    dim3 dims = make_uint3(nx-1, ny, nz);
    pass1gpu<<<dims, 1>>>(pointValues, isoval, edgeCases);
}
*/
///////////////////////////////////////////////////////////////////////////////
// *Time: 8.2e-5
///////////////////////////////////////////////////////////////////////////////

/*
__global__
void pass1gpu1(
    scalar_t* pointValues,  // input
    scalar_t isoval,        // input
    uchar* edgeCases)       // output
{
    int i = blockIdx.x;
    int j = blockIdx.y;
    int k = blockIdx.z;

    int ny = blockDim.y;
    int nx = blockDim.x + 1;

    edgeCases[k*ny*(nx-1) + j*(nx-1) + i] =
        calcCaseEdge(
            pointValues[k*ny*nx + j*nx + i + 0] >= isoval,
            pointValues[k*ny*nx + j*nx + i + 1] >= isoval);
}

__global__
void pass1gpu2(
    uchar* edgeCases,    // input
    int nx,              // input
    FlyingEdgesAlgorithm::gridEdge* gridEdges) // output
{
    int j = blockIdx.y;
    int k = blockIdx.z;

    int ny = blockDim.y;

    FlyingEdgesAlgorithm::gridEdge& grid = gridEdges[k*ny + j];

    for(int i = 0; i != nx-1; ++i)
    {
        uchar const& edge = edgeCases[k*ny*(nx-1) + j*(nx-1) + i];
        if(edge == 1 || edge == 2)
        {
            grid.xl = i;
            break;
        }
    }

    for(int i = nx-2; i != -1; ++i)
    {
        uchar const& edge = edgeCases[k*ny*(nx-1) + j*(nx-1) + i];
        if(edge == 1 || edge == 2)
        {
            grid.xr = i;
            break;
        }
    }
}

__global__
void pass1gpu222(
    scalar_t* pointValues,                     // input
    scalar_t isoval,                           // input
    int nx,                                    // input
    int ny,
    uchar* edgeCases,                          // output
    FlyingEdgesAlgorithm::gridEdge* gridEdges) // output
{
//    int j = blockIdx.y;
//    int k = blockIdx.z;
//
//    int ny = blockDim.y;
    int j = threadIdx.y;
    int k = threadIdx.z;

    scalar_t* curPointValues = pointValues + k*nx*ny + j*nx;
    uchar* curEdgeCases = edgeCases + k*(nx-1)*ny + j*(nx-1);
    FlyingEdgesAlgorithm::gridEdge& curGridEdge = gridEdges[k*ny + j];

    bool isGE[2];
    isGE[0] = (curPointValues[0] >= isoval);
    for(int i = 1; i != nx; ++i)
    {
        isGE[i%2] = (curPointValues[i] >= isoval);
        curEdgeCases[i-1] = calcCaseEdge(isGE[(i+1)%2], isGE[i%2]);

        if(curEdgeCases[i-1] == 1 || curEdgeCases[i-1] == 2)
        {
            if(curGridEdge.xl == 0)
                curGridEdge.xl = i-1;
            curGridEdge.xr = i;
        }
    }
}

__global__
void pass1gpu333(
    scalar_t* pointValues,                     // input
    scalar_t isoval,                           // input
    int nx, int ny, int nz,                    // input
    uchar* edgeCases,                          // output
    FlyingEdgesAlgorithm::gridEdge* gridEdges) // output
{
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (j >= ny || k >= nz)
        return;

    scalar_t* curPointValues = pointValues + k*nx*ny + j*nx;
    uchar* curEdgeCases = edgeCases + k*(nx-1)*ny + j*(nx-1);
    FlyingEdgesAlgorithm::gridEdge& curGridEdge = gridEdges[k*ny + j];

    bool isGE[2];
    isGE[0] = (curPointValues[0] >= isoval);
//    for(int i = 1; i != nx; ++i)
//    {
//        isGE[i%2] = (curPointValues[i] >= isoval);
//        curEdgeCases[i-1] = calcCaseEdge(isGE[(i+1)%2], isGE[i%2]);
//
//        if(curEdgeCases[i-1] == 1 || curEdgeCases[i-1] == 2)
//        {
//            if(curGridEdge.xl == 0)
//                curGridEdge.xl = i-1;
//            curGridEdge.xr = i;
//        }
//    }
}
*/
