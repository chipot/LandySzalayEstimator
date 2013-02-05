#include "log.hh"
#include <iostream>

namespace llog {

notifier<level_notice> notice(std::cout);
notifier<level_warning> warn(std::clog);
notifier<level_fatal> fatal(std::clog);
notifier<level_debug> debug(std::cerr);
    
} /* llog */

