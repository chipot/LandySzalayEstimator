#pragma once

namespace htm {

struct trixel;

struct PointInfo
{
    double          _ra;       // Right ascension
    double          _dec;      // Declination
    struct trixel  *_current;  // Current position in HTM
};

} // htm
