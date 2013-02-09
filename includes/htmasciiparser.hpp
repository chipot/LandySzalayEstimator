#pragma once

#include <tbb/concurrent_vector.h>
#include <utility>

namespace htm {

class HTM;

class HTMAsciiParser
{
 public:
  tbb::concurrent_vector<std::pair<double, double>> Parse(std::string& filename);

  HTMAsciiParser() = default;
  ~HTMAsciiParser() = default;
};

}
