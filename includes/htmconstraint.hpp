#ifndef __ICODF_HTM_HTMCONSTRAINT__H__
#define __ICODF_HTM_HTMCONSTRAINT__H__

// C++ includes
#include <vector>

namespace htm {

struct trixel;

struct Constraint
{
    std::vector<trixel*> _inside;
    std::vector<trixel*> _partial;
};

} // ICoDF_HTM

#endif /* __ICODF_HTM_HTMCONSTRAINT__H__ */
