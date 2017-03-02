#ifndef FLYINGEDGES_CONFIG_H_
#define FLYINGEDGES_CONFIG_H_

#include <array>

using std::size_t;

using scalar_t = float;

using uchar = unsigned char;

using cube_t = std::array<std::array<scalar_t, 3>, 8>;
using scalarCube_t = std::array<scalar_t, 8>;

#endif
