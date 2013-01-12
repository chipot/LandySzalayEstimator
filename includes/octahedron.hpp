#pragma once

#include <vector>
#include "trixel.hpp"

namespace ICoDF_HTM {

typedef struct Octahedron_s
{
    std::string _name;
    trixel_t** _rootTrixels; // [0] = S0, [4] = N0
} Octahedron_t;

} // ICoDF_HTM
