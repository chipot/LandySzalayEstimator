#pragma once

namespace ICoDF_HTM {

struct trixel_s;

typedef struct PointInfo_s
{
    double _ra;                 // Right ascension
    double _dec;                // Declination
    struct trixel_s* _current;  // Current position in HTM
} PointInfo_t;

} // ICoDF_HTM
