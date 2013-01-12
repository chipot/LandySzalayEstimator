#ifndef __BLINK_HTM_TRIXEL_HPP__
#define __BLINK_HTM_TRIXEL_HPP__

// C Includes
#include <assert.h>

// C++ Includes
#include <string>
#include <sstream>

// Eigen Includes
#include <Eigen/Dense>

// ICoDF
#include "logservice.hpp"
#include "pointinfo.hpp"

using namespace ICoDF;
using namespace ICoDF_HTM;

namespace ICoDF_HTM {

struct PointInfo_s;

/// Structure that define a trixel (htm base object)
typedef struct trixel_s
{
    struct trixel_s**       _children;  //< trixel's subtrixels
    Eigen::Vector3d*        _vertices;  //< Trixel's vertices
    bool                    _reverse;   //< is an upside-down trixel ?
    std::string             _HTMId;     //< N10120112121101
    unsigned int            _nbChildObject; //< Number of objects contained in this trixel.
    struct PointInfo_s*     _info;      //< Point information structure for the actual 
} trixel_t;

/// Delte the given trixel (but not the PointInfo_t)
void ClearTrixel(trixel_t* trixel);

/// Delete a trixel's children
void ClearTrixelChildren(trixel_t *parent);

/// Compute a trixel's midpoint vectors
Eigen::Vector3d* ComputeTrixelMidpoints(trixel_t* trixel);

/// Create A Root Trixel
trixel_t* CreateRootTrixel(std::string HTMId);

/// Create a new trixel from its parent information and index
trixel_t* CreateTrixelChild(trixel_t *parent, unsigned short int& index);

/// Create a container for the 4 trixels that correspond to a new level 
trixel_t** CreateTrixelChildren(trixel_t* parent);

/// Return which subtrixel the point is (ra dec version)
unsigned short int GetIndex(trixel_t* trixel, double& ra, double& dec);

/// Return which subtrixel the point is (vector version)
unsigned short int GetIndex(trixel_t* trixel, Eigen::Vector3d& pointVector);

/// Return which subtrixel the point is (object version)
unsigned short int GetIndex(trixel_t* trixel, PointInfo_t* object);

/// Init a new trixel with default values
void InitTrixel(trixel_t* trixel);

/// Check if the given right ascension value is correct
bool IsCorrectRA(double& ra);

/// Check if the given declination value is correct
bool IsCorrectDEC(double& dec);

} // ICoDF_HTM

#endif /* !__BLINK_HTM_TRIXEL_HPP__ */
