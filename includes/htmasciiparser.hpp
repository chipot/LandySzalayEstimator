#pragma once

#include <string>
#include <vector>
#include <utility>

namespace htm {

class HTM;

class HTMAsciiParser
{
 public:
  void Parse(std::string& filename);
  void Populate();
  void UniformNumberGenerator(const double& raMin, const double& raMax, const double& decMin, const double& decMax);
  unsigned int& getNbObj(void);
  // DEFAULT CTOR
  HTMAsciiParser(HTM *);

  // DEFAULT DTOR
  ~HTMAsciiParser();
 private:
  HTM* _htm;
  std::vector<std::pair<double, double>> _cache;
  unsigned int nbObj;
};

}
