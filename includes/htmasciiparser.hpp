#pragma once

#include <vector>
#include <utility>

namespace htm {

class HTM;

class HTMAsciiParser
{
 public:
  std::vector<std::pair<double, double>> Parse(std::string& filename);

  HTMAsciiParser() = default;
  ~HTMAsciiParser() = default;
};

}
