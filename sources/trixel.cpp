#include "../includes/trixel.hpp"

namespace /* annon */ {

template <typename T>
constexpr T
max()
{
    return (T)~0;
} /* annon */

}

// ------------------------------------------------------------------------------
// CREATETRIXELCHILDREN
trixel_t** ICoDF_HTM::CreateTrixelChildren(trixel_t *parent)
{
    if (parent->_children != NULL)
    {
        for (int i = 0; i < 4; ++i)
        {
            if (parent->_children[i] != NULL)
            {
                LS_ADDMSG(LogService::NOTICE, "ICoDF::CreateTrixelChildren", "Trixel already have child(ren)");
            }
        }
    }
    else
    {
        parent->_children = new trixel_t*[4];
        for (int i = 0; i <  4; ++i)
        {
            parent->_children[i] = NULL;
        }
    }
    return parent->_children;
}

// -------------------------------------------------------------------------------
Eigen::Vector3d* ICoDF_HTM::ComputeTrixelMidpoints(trixel_t* trixel)
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
trixel_t* ICoDF_HTM::CreateRootTrixel(std::string HTMId)
{
    trixel_t* trixel = new trixel_t();
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
trixel_t* ICoDF_HTM::CreateTrixelChild(trixel_t* parent, unsigned short int& index)
{
    if (parent->_children == NULL)
    {
        LS_ADDMSG(LogService::NOTICE, "ICoDF::CreateTrixelChild", "Trixel as no container for children");
        CreateTrixelChildren(parent);
    }

    assert(index < 4);

    if (parent->_children[index] == NULL)
    {
        static std::stringstream tmp;
        parent->_children[index] = new trixel_t();
        InitTrixel(parent->_children[index]);
        tmp.str("");
        tmp << parent->_HTMId << index;
        parent->_children[index]->_HTMId = tmp.str();
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
                LS_ADDMSG(LogService::FATAL, "ICoDF::CreateTrixelChild", "Given <index> is out of bound");
                delete [] midPoints; 
                return NULL;
        }
        delete [] midPoints; 
    }
    else
    {
        std::stringstream tmp;
        tmp << "SubTrixel [" << parent->_HTMId << index << "] already exists";
        LS_ADDMSG(LogService::NOTICE, "ICoDF::CreateTrixelChild", tmp.str());
    }

    return parent->_children[index];
}

// -------------------------------------------------------------------
// CLEARTRIXELCHILDREN
void ICoDF_HTM::ClearTrixelChildren(trixel_t *parent)
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
void ICoDF_HTM::ClearTrixel(trixel_t *trixel)
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
unsigned short int ICoDF_HTM::GetIndex(trixel_t* trixel, double& ra, double& dec)
{
    if (IsCorrectRA(ra) && IsCorrectDEC(dec))
    {
        double rProjection = sin(90 - abs(dec));
        double x = rProjection * cos(ra);
        double y = rProjection * sin(ra);
        double z = cos(90 - abs(dec));
        Eigen::Vector3d p(x, y, z);
        unsigned short int ret = GetIndex(trixel, p);
        if (ret == max<unsigned short>())
        {
            std::stringstream tmp;
            tmp << " T0 " << std::endl << trixel->_vertices[0] << std::endl
                << " T1 " << std::endl << trixel->_vertices[1] << std::endl
                << " T2 " << std::endl << trixel->_vertices[2] << std::endl
                << " P " << std::endl << p << std::endl;
        }
        return ret;
    }
    else
    {
        std::stringstream tmp;
        tmp << "Given <ra> [" << ra << "] or <dec> [" << dec << "] is out of bounds";
        LS_ADDMSG(LogService::WARNING, "ICoDF_HTM", tmp.str());
    }
    return max<unsigned short>();
}

// --------------------------------------------------------------------
// GETINDEX (vector version)
unsigned short int ICoDF_HTM::GetIndex(trixel_t* trixel, Eigen::Vector3d& p)
{
    if (trixel != NULL && NULL != trixel->_vertices)
    {
        unsigned short int index = max<unsigned short int>();

        Eigen::Vector3d* v = trixel->_vertices;
        Eigen::Vector3d* w = ComputeTrixelMidpoints(trixel);

        // HERE
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
            std::cout << "Incorrect : " << trixel->_HTMId << std::endl << "--- v1" << std::endl <<  trixel->_vertices[0] << std::endl << "--- v2" << std::endl << trixel->_vertices[1] << std::endl << "--- v3" << std::endl << trixel->_vertices[2] << std::endl << "--- p" << std::endl << p << std::endl;

        delete [] w;
        return index;
    }
    else
    {
        LS_ADDMSG(LogService::WARNING, "ICoDF_HTM::GetIndex", "Given <trixel> or its vertices has a NULL value");
        return max<unsigned short>();
    }
}

// --------------------------------------------------------------------
// GETINDEX (PointInfo version)
unsigned short int ICoDF_HTM::GetIndex(trixel_t* trixel, PointInfo_t* pointInfo)
{
    unsigned short int index = GetIndex(trixel, pointInfo->_ra, pointInfo->_dec);
    return index;
}

// --------------------------------------------------------------------
// INITTRIXEL
void ICoDF_HTM::InitTrixel(trixel_t* trixel)
{
    if (trixel == NULL)
    {
        LS_ADDMSG(LogService::WARNING, "ICoDF_HTM::InitTrixel", "Given <trixel> has a NULL value");
    }
    else
    {
        trixel->_children = NULL;
        trixel->_vertices = NULL;
        trixel->_vertices = new Eigen::Vector3d[3];
        trixel->_reverse = false;
        trixel->_HTMId = "";
        trixel->_nbChildObject = 0;
        trixel->_info = NULL;      
    }
}

// --------------------------------------------------------------------
// ISCORRECTRA
bool ICoDF_HTM::IsCorrectRA(double& ra)
{
    if (ra >= 0 && ra < 360)
        return true;
    return false;
}

// --------------------------------------------------------------------
// ISCORRECTDEC
bool ICoDF_HTM::IsCorrectDEC(double& dec)
{
    if (dec > -90 && dec <= 90)
        return true;
    return false;
}
