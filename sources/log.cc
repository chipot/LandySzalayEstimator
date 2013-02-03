#include "log.hh"
#include <iostream>

namespace llog {

notifier<level_notice> notice(std::clog);
notifier<level_warning> warn(std::clog);
notifier<level_fatal> fatal(std::clog);
    
} /* llog */

