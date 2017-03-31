#ifndef THRUST_CONFIG_H
#define THRUST_CONFIG_H

using scalar_t = float;
using uchar = unsigned char;

template <typename T>
struct t3
{
//    t3()
//    {}
//
//    t3(t3 const& other)
//      : x(x), y(y), z(z)
//    {}
//
    t3(T const& x, T const& y, T const& z)
      : x(x), y(y), z(z)
    {}

    T x;
    T y;
    T z;
};

#endif
