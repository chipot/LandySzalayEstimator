#ifndef __BLINK_HTM_TRIXEL_HPP__
#define __BLINK_HTM_TRIXEL_HPP__

#include <string>
#include <Eigen/Dense>

namespace htm {

struct PointInfo;

/// Structure that define a trixel (htm base object)
struct trixel
{
    unsigned int            _nbChildObject; //< Number of objects contained in this trixel.
    struct trixel         **_children;  //< trixel's subtrixels
    Eigen::Vector3d        *_vertices;  //< Trixel's vertices
    struct PointInfo       *_info;      //< Point information structure for the actual 
    std::string             _HTMId;     //< N10120112121101
};

/// Delte the given trixel (but not the PointInfo)
void ClearTrixel(trixel* trixel);

/// Delete a trixel's children
void ClearTrixelChildren(trixel *parent);

/// Compute a trixel's midpoint vectors
Eigen::Vector3d* ComputeTrixelMidpoints(trixel* trixel);

/// Create A Root Trixel
trixel* CreateRootTrixel(std::string HTMId);

/// Create a new trixel from its parent information and index
trixel* CreateTrixelChild(trixel *parent, unsigned short int& index);

/// Create a container for the 4 trixels that correspond to a new level 
trixel** CreateTrixelChildren(trixel* parent);

/// Return which subtrixel the point is (ra dec version)
unsigned short int GetIndex(trixel* trixel, double& ra, double& dec);

/// Return which subtrixel the point is (vector version)
unsigned short int GetIndex(trixel* trixel, Eigen::Vector3d& pointVector);

/// Return which subtrixel the point is (object version)
unsigned short int GetIndex(trixel* trixel, PointInfo* object);

/// Init a new trixel with default values
void InitTrixel(trixel* trixel);

/// Check if the given right ascension value is correct
bool IsCorrectRA(double const& ra);

/// Check if the given declination value is correct
bool IsCorrectDEC(double const& dec);

} // htm

#endif /* !__BLINK_HTM_TRIXEL_HPP__ */
