#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include <thrust/functional.h>

#include <thrust/sequence.h>
#include <thrust/transform.h>
#include <thrust/for_each.h>
#include <thrust/copy.h>

#include <thrust/tuple.h>

#include <thrust/iterator/zip_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>

#include <iostream>
#include <string.h>
#include <cstdlib>

#include "../util/LoadImage.h"
#include "../util/SaveTriangleMesh.h"

#include "../util/Timer.h"
#include "../mantevoCommon/YAML_Doc.hpp"

using scalar_t = float;
using uchar = unsigned char;

/*
 *
 * Prints out 0, 5, 10, ..., 95
 *
 *
    thrust::device_vector<int> dv(100);
    thrust::sequence(dv.begin(), dv.end());

    access_iter_helper<int> aih(dv.begin(), 20, 5);

    for(typename access_iter_helper<int>::iterator i = aih.begin(); i != aih.end(); ++i)
    {
        std::cout << (*i)[0] << std::endl;
    }
*/
template <typename T>
struct iter_access_helper
{
    struct iter_access_op
      : public thrust::unary_function<
            int,
            typename thrust::device_vector<T>::iterator>
    {
        iter_access_op(
            typename thrust::device_vector<T>::iterator const begin,
                int const& n)
          : begin(begin),
            n(n)
            {}

        __host__ __device__
        typename thrust::device_vector<T>::iterator
        operator()(int x) const
        {
            return begin + n*x;
        }

        typename thrust::device_vector<T>::iterator const begin;
        int const n;
    };

    using counting_iterator = thrust::counting_iterator<int>;
    using iterator =
        thrust::transform_iterator<
            iter_access_op,
            counting_iterator>;

    iter_access_helper(
        typename thrust::device_vector<T>::iterator const begin,
        int numTotal,
        int step)
      : begin_iter(
            iterator(
                thrust::make_counting_iterator(0),
                iter_access_op(begin, step))),
        end_iter(begin_iter + numTotal)
    {}

    iterator begin() const
    {
        return begin_iter;
    }

    iterator end() const
    {
        return end_iter;
    }

private:
    iterator begin_iter;
    iterator end_iter;
};

template <typename T>
struct t3
{
    t3()
    {}

    t3(T const& x, T const& y, T const& z)
      : x(x), y(y), z(z)
    {}

    T x;
    T y;
    T z;
};

struct edge_case
  : public thrust::unary_function<
        thrust::tuple<bool, bool>,
        uchar>
{
    __host__ __device__
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


struct grid_change
  : public thrust::unary_function<
        int,
        int>
{
    grid_change(t3<const int> n)
      : n(n)
    {}

    __host__ __device__
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

struct isGE
  : public thrust::unary_function<
        scalar_t,
        bool>
{
    isGE(scalar_t const& isoval)
      : isoval(isoval)
    {}

    __host__ __device__
    bool
    operator()(scalar_t const& val)
    {
        return val >= isoval;
    }

    scalar_t const isoval;
};

struct useless
  : public thrust::unary_function<
        thrust::tuple<bool, bool>,
        void>
{
    __host__ __device__
    void
    operator()(thrust::tuple<bool, bool> const& b) const
    {
    }
};

struct multiply_transform
  : public thrust::unary_function<int, int>
{
    multiply_transform(int n)
      : n(n)
    {}

    multiply_transform()
      : n(1)
    {}

    __host__ __device__
    int operator()(int x) const
    {
        return n*x;
    }
private:
    int n;
};

struct useless8
  : public thrust::unary_function<
        scalar_t,
        uchar>
{
    __host__ __device__
    uchar operator()(scalar_t const& s) const
    {
        return 1;
    }
};

struct useless9
  : public thrust::unary_function<
        thrust::tuple<scalar_t, uchar>,
        void>
{
    __host__ __device__
    void operator()(thrust::tuple<scalar_t, uchar> const& vvv) const
    {

    }
};

struct useless10
  : public thrust::unary_function<
        scalar_t,
        void>
{
    __host__ __device__
    void operator()(scalar_t const& w)
    {
    }
};

struct Image3D
{
public:
    Image3D(thrust::host_vector<scalar_t> data_,
            scalar_t const& spacingX,
            scalar_t const& spacingY,
            scalar_t const& spacingZ,
            scalar_t const& zeroPosX,
            scalar_t const& zeroPosY,
            scalar_t const& zeroPosZ,
            int const& nx,
            int const& ny,
            int const& nz)
      : data(data_),
        spacing(spacingX, spacingY, spacingZ),
        zeroPos(zeroPosX, zeroPosY, zeroPosZ),
        n(nx, ny, nz),
        ray_helper(data.begin(), n.y*n.z, n.x)
    {}

    using ray_iterator =
        typename iter_access_helper<scalar_t>::iterator;

    ray_iterator ray_begin() const
    {
        return ray_helper.begin();
    }

    ray_iterator ray_end() const
    {
        return ray_helper.end();
    }

    int xdimension() const { return n.x; }
    int ydimension() const { return n.y; }
    int zdimension() const { return n.z; }

    using iterator =
        typename thrust::device_vector<scalar_t>::const_iterator;

    iterator begin() const
    {
        return data.begin();
    }

    iterator end() const
    {
        return data.end();
    }

private:
    thrust::device_vector<scalar_t> data;
public:
    t3<const scalar_t> spacing;
    t3<const scalar_t> zeroPos;
    t3<const int> n;
private:
    iter_access_helper<scalar_t> ray_helper;
};


struct FlyingEdges
{
    FlyingEdges(Image3D const& image, scalar_t const& isoval)
      : image(image),
        isoval(isoval),
        n(image.n),
        left(n.y*n.z, n.x),
        right(n.y*n.z, 0),
        edgeCases((n.x-1)*n.y*n.z)
    {}

    void pass1()
    {
        // set values of trim.edgeCases
/*
        using isGE_iter_t =
            thrust::transform_iterator<
                isGE,
                typename Image3D::iterator>;

        isGE_iter_t isGE_iter =
            thrust::make_transform_iterator(
                image.begin(),
                isGE(isoval));

        using element_iterator_t =
            thrust::zip_iterator<
                thrust::tuple<
                    isGE_iter_t,
                    isGE_iter_t> >;

        element_iterator_t element_iterator =
            thrust::make_zip_iterator(
                thrust::make_tuple(
                    isGE_iter,
                    isGE_iter + 1));

        using index_iterator_t =
            thrust::transform_iterator<
                grid_change,
                thrust::counting_iterator<int> >;

        index_iterator_t index_iterator =
            thrust::make_transform_iterator(
                thrust::make_counting_iterator(0),
                grid_change(n));

        using permutation_iterator_t =
            thrust::permutation_iterator<
                element_iterator_t,
                index_iterator_t>;

        permutation_iterator_t input_iter_begin =
            thrust::make_permutation_iterator(
                element_iterator,
                index_iterator);

        permutation_iterator_t input_iter_end =
            input_iter_begin + (n.x-1)*n.y*n.z;
*/
//        thrust::for_each(
//            input_iter_begin,
//            input_iter_end,
//            useless());

        auto beg_wtf = thrust::make_zip_iterator(
            thrust::make_tuple(
                image.begin(),
                edgeCases.begin()));

        auto end_wtf = beg_wtf + 100;

//        thrust::for_each(
//            beg_wtf,
//            end_wtf,
//            useless9());

//        thrust::transform(
//            image.begin(),
//            image.end(),
//            edgeCases.begin(),
//            useless8());

//        thrust::transform(
//            input_iter_begin,
//            input_iter_end,
//            edgeCases.begin(),
//            edge_case());

        int count = 0;
        int count0 = 0;
        int count1 = 0;
        int count2 = 0;
        int count3 = 0;
        for(uchar const& c: edgeCases)
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

        // set values of xl and xr

        // TODO TODO TODO TODO TODO TODO
        // typedef host vector and iterator and hsit
        // what happens with big vector?
        // First:  organize code, move most of this code to old file
        // Second: come up with way to set xl and xr.

    }
private:
    Image3D const& image;
    scalar_t const isoval;
    t3<const int> n;

    thrust::host_vector<int> left;
    thrust::host_vector<int> right;
    thrust::host_vector<uchar> edgeCases;
};

int main(int argc, char* argv[])
{
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

    std::vector<scalar_t> data;
    scalar_t spacingX;
    scalar_t spacingY;
    scalar_t spacingZ;
    scalar_t zeroPosX;
    scalar_t zeroPosY;
    scalar_t zeroPosZ;
    int nx;
    int ny;
    int nz;

    util::loadImage_thrust(
        vtkFile,
        data,
        spacingX,
        spacingY,
        spacingZ,
        zeroPosX,
        zeroPosY,
        zeroPosZ,
        nx,
        ny,
        nz);

    int nx_actually = 100;
    int sz = 512*512*nx_actually;
    std::vector<scalar_t> cut_data(sz);
    std::copy(data.begin(), data.begin() + sz, cut_data.begin());

    thrust::host_vector<scalar_t> host_data(cut_data);

    Image3D image(
        host_data,
        spacingX,
        spacingY,
        spacingZ,
        zeroPosX,
        zeroPosY,
        zeroPosZ,
        nx_actually,
        ny,
        nz);

    doc.add("File x-dimension", image.xdimension());
    doc.add("File y-dimension", image.ydimension());
    doc.add("File z-dimension", image.zdimension());

    util::Timer runTime;

    util::Timer toDeviceTime;
    FlyingEdges algo(image, isoval);
    toDeviceTime.stop();

    util::Timer runTimePass1;
    algo.pass1();
    runTimePass1.stop();

    util::Timer runTimePass2;
//    algo.pass2();
    runTimePass2.stop();

    util::Timer runTimePass3;
//    algo.pass3();
    runTimePass3.stop();

    util::Timer runTimePass4;
//    algo.pass4();
    runTimePass4.stop();

    util::Timer fromDeviceTime;
//    util::TriangleMesh mesh = algo.moveOutput();
    fromDeviceTime.stop();

    runTime.stop();

//    doc.add("Number of vertices in mesh", mesh.numberOfVertices());
//    doc.add("Number of triangles in mesh", mesh.numberOfTriangles());

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
/*
    thrust::device_vector<int> dv(100);
    thrust::sequence(dv.begin(), dv.end());

    using citer = thrust::counting_iterator<int>;
    using titer =
        thrust::transform_iterator<
            iter_access<int>,
            citer>;

    titer t_beg = titer(
        thrust::make_counting_iterator(0),
        iter_access<int>(dv.begin(), 10));

    titer t_end = t_beg + 10;

    for(; t_beg != t_end; ++t_beg)
    {
        std::cout << (*t_beg)[0] << std::endl;
    }
*/
}

