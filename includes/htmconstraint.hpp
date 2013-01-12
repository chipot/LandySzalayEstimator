#ifndef __ICODF_HTM_HTMCONSTRAINT__H__
#define __ICODF_HTM_HTMCONSTRAINT__H__

// C++ includes
#include <list>

// ICoDF includes
#include "trixel.hpp"

namespace ICoDF_HTM {

typedef struct HTMConstraint_s
{
    std::list<trixel_t*> _inside;
    std::list<trixel_t*> _partial;
} HTMConstraint_t;

} // ICoDF_HTM

#endif /* __ICODF_HTM_HTMCONSTRAINT__H__ */
