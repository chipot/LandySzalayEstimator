#pragma once

#include <vector>
#include <utility>

namespace htm {

class HTM;

class HTMAsciiParser
{
 public:
  std::vector<std::pair<double, double>> Parse(std::string& filename);

  // DEFAULT CTOR
  HTMAsciiParser(HTM *);

  // DEFAULT DTOR
  ~HTMAsciiParser();
 private:
  HTM* _htm;
};

}
