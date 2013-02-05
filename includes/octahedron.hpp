#pragma once
#ifndef OCTAHEDRON_6769ZRSD
#define OCTAHEDRON_6769ZRSD

namespace htm {

struct trixel;

struct Octahedron
{
    std::string _name;
    trixel** _rootTrixels; 
    // [0] = S0, [4] = N0
    // [1] = S1, [5] = N1
    // [2] = S2, [6] = N2
    // [3] = S3, [7] = N3
};

}

#endif /* end of include guard: OCTAHEDRON_6769ZRSD */

