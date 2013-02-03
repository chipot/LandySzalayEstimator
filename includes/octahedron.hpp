#pragma once
#ifndef OCTAHEDRON_6769ZRSD
#define OCTAHEDRON_6769ZRSD

namespace htm {

struct trixel;

struct Octahedron
{
    std::string _name;
    trixel** _rootTrixels; // [0] = S0, [4] = N0
};

}

#endif /* end of include guard: OCTAHEDRON_6769ZRSD */

