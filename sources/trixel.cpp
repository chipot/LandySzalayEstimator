// C Includes
#include <assert.h>

// C++ Includes
#include <iostream>
#include <sstream>

// Eigen Includes
#include <Eigen/Dense>

// ICoDF
#include "log.hh"
#include "pointinfo.hpp"
#include "trixel.hpp"

namespace /* annon */ {

template <typename T>
constexpr T
max()
{
    return (T)~0;
} /* annon */

}

namespace htm {

// ------------------------------------------------------------------------------
// CREATETRIXELCHILDREN
trixel** CreateTrixelChildren(trixel *parent)
{
    if (parent->_children != NULL)
    {
        for (int i = 0; i < 4; ++i)
        {
            if (parent->_children[i] != NULL)
            {
                llog::debug["ICoDF::CreateTrixelChildren"]
                    <<  "Trixel already have child(ren)" << std::endl;
            }
        }
    }
    else
    {
        parent->_children = new trixel*[4];
        for (int i = 0; i <  4; ++i)
        {
            parent->_children[i] = NULL;
        }
    }
    return parent->_children;
}

// -------------------------------------------------------------------------------
Eigen::Vector3d* ComputeTrixelMidpoints(trixel* trixel)
{
    static Eigen::Vector3d tmp;
    Eigen::Vector3d* midPoints = new Eigen::Vector3d[3];

    tmp = trixel->_vertices[1] + trixel->_vertices[2];
    midPoints[0] = tmp / tmp.norm();
    tmp = trixel->_vertices[0] + trixel->_vertices[2];
    midPoints[1] = tmp / tmp.norm();
    tmp = trixel->_vertices[0] + trixel->_vertices[1];
    midPoints[2] = tmp / tmp.norm();

    return midPoints;
}

// -------------------------------------------------------------------------------
// CREATEROOTTRIXEL
trixel* CreateRootTrixel(std::string HTMId)
{
    trixel* trixel = new struct trixel();
    InitTrixel(trixel);
    trixel->_HTMId = HTMId;
    return trixel;
}

// -------------------------------------------------------------------------------
// CREATETRIXELCHILD
// Check if trixel children structure already exists
// Assert for given <index>
// Check if select child does not exist
//   Create new vertices on side midpoints
//   Set vertices that defines the new subtrixel
// else display a message
// return the pointer
trixel* CreateTrixelChild(trixel* parent, unsigned short int& index)
{
    if (parent->_children == NULL)
    {
        llog::debug["ICoDF::CreateTrixelChild"]
            <<  "Trixel has no container for children" << std::endl;
        CreateTrixelChildren(parent);
    }

    assert(index < 4);

    if (parent->_children[index] == NULL)
    {
        parent->_children[index] = new trixel();
        InitTrixel(parent->_children[index]);
        parent->_children[index]->_HTMId = parent->_HTMId + std::to_string(index);
        Eigen::Vector3d* midPoints = ComputeTrixelMidpoints(parent);

        switch (index)
        {
            case 0:
                parent->_children[index]->_vertices[0] = parent->_vertices[0];
                parent->_children[index]->_vertices[1] = midPoints[2];
                parent->_children[index]->_vertices[2] = midPoints[1];
                break;
            case 1:
                parent->_children[index]->_vertices[0] = parent->_vertices[1];
                parent->_children[index]->_vertices[1] = midPoints[0];
                parent->_children[index]->_vertices[2] = midPoints[2];
                break;
            case 2:
                parent->_children[index]->_vertices[0] = parent->_vertices[2];
                parent->_children[index]->_vertices[1] = midPoints[1];
                parent->_children[index]->_vertices[2] = midPoints[0];
                break;
            case 3:
                parent->_children[index]->_vertices[0] = midPoints[0];
                parent->_children[index]->_vertices[1] = midPoints[1];
                parent->_children[index]->_vertices[2] = midPoints[2];
                break;
            default:
                llog::warn["ICoDF::CreateTrixelChild"]
                    << "Given <index> is out of bound" << std::endl;
                delete [] midPoints; 
                return NULL;
        }
        delete [] midPoints; 
    }
    else
    {
        llog::warn["ICoDF::CreateTrixelChild"]
            << "SubTrixel [" << parent->_HTMId << index << "] already exists"
            << std::endl;
    }

    return parent->_children[index];
}

// -------------------------------------------------------------------
// CLEARTRIXELCHILDREN
void ClearTrixelChildren(trixel *parent)
{
    if (parent->_children != NULL)
    {
        for (short int i = 0; i < 4; ++i)
        {
            if (parent->_children[i] != NULL)
            {
                ClearTrixel(parent->_children[i]);
                parent->_children[i] = NULL;
            }
        }
        delete [] parent->_children;
        parent->_children = NULL;
    }
}

// ---------------------------------------------------------------------
// CLEARTRIXELCHILDREN
void ClearTrixel(trixel *trixel)
{
    if (trixel->_vertices != NULL)
    {
        delete [] trixel->_vertices;
        trixel->_vertices = NULL;
    }
    ClearTrixelChildren(trixel);
}

// ---------------------------------------------------------------------
// GETINDEX (astral coordinates version)
unsigned short int GetIndex(trixel* trixel, double& ra, double& dec)
{
    if (IsCorrectRA(ra) && IsCorrectDEC(dec))
    {
        double rProjection = sin(90 - abs(dec));
        double x = rProjection * cos(ra);
        double y = rProjection * sin(ra);
        double z = cos(90 - abs(dec));
        Eigen::Vector3d p(x, y, z);
        unsigned short int ret = GetIndex(trixel, p);
        if (ret == max<unsigned short>()) // Is it not good ..?
        {
            // Did nothing in the first place...
        }
        return ret;
    }
    else
    {
        llog::warn["HTM"]
            << "Given <ra> [" << ra << "] or <dec> [" << dec << "] is out of bounds"
            << std::endl;
    }
    return max<unsigned short>();
}

// --------------------------------------------------------------------
// GETINDEX (vector version)
unsigned short int GetIndex(trixel* trixel, Eigen::Vector3d& p)
{
    if (trixel != NULL && NULL != trixel->_vertices)
    {
        unsigned short int index = max<unsigned short int>();

        Eigen::Vector3d* v = trixel->_vertices;
        Eigen::Vector3d* w = ComputeTrixelMidpoints(trixel);

        // HERE // Here WHAT ??!
        if (v[0].cross(w[2]).dot(p) > 0 &&
            w[2].cross(w[1]).dot(p) > 0 &&
            w[1].cross(v[0]).dot(p) > 0)
            index = 0;
        else if (v[1].cross(w[0]).dot(p) > 0 &&
                 w[0].cross(w[2]).dot(p) > 0 &&
                 w[2].cross(v[1]).dot(p) > 0)
            index = 1;
        else if (v[2].cross(w[1]).dot(p) > 0 &&
                 w[1].cross(w[0]).dot(p) > 0 &&
                 w[0].cross(v[2]).dot(p) > 0)
            index = 2;
        else if (w[0].cross(w[1]).dot(p) > 0 &&
                 w[1].cross(w[2]).dot(p) > 0 &&
                 w[2].cross(w[0]).dot(p) > 0)
            index = 3;

        if (index == max<unsigned short>())
            // Shit is incorrect, but why, what, how ? God only knows
            llog::warn["GetIndex"]  << "Incorrect : " << trixel->_HTMId << std::endl
                << "--- v1" << std::endl <<  trixel->_vertices[0] << std::endl
                << "--- v2" << std::endl << trixel->_vertices[1] << std::endl
                << "--- v3" << std::endl << trixel->_vertices[2] << std::endl
                << "--- p" << std::endl << p << std::endl;

        delete [] w;
        return index;
    }
    else
    {
        llog::warn["GetIndex"]
            <<  "Given <trixel> or its vertices has a NULL value" << std::endl;
        return max<unsigned short>();
    }
}

// --------------------------------------------------------------------
// GETINDEX (PointInfo version)
unsigned short int GetIndex(trixel* trixel, PointInfo* pointInfo)
{
    unsigned short int index = GetIndex(trixel, pointInfo->_ra, pointInfo->_dec);
    return index;
}

// --------------------------------------------------------------------
// INITTRIXEL
void InitTrixel(trixel* trixel)
{
    if (trixel == NULL)
    {
        llog::warn["InitTrixel"]
            <<  "Given <trixel> has a NULL value" << std::endl;
    }
    else
    {
        trixel->_children = NULL;
        trixel->_vertices = NULL;
        trixel->_vertices = new Eigen::Vector3d[3];
        trixel->_HTMId = "";
        trixel->_nbChildObject = 0;
        trixel->_info = NULL;      
    }
}

// --------------------------------------------------------------------
// ISCORRECTRA
bool IsCorrectRA(double const& ra)
{
    if (ra >= 0 && ra < 360)
        return true;
    return false;
}

// --------------------------------------------------------------------
// ISCORRECTDEC
bool IsCorrectDEC(double const& dec)
{
    if (dec > -90 && dec <= 90)
        return true;
    return false;
}

}
