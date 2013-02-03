#pragma once

#include <string>

namespace htm {

class HTM;

class HTMAsciiParser
{
 public:
  void Parse(std::string& filename);
  void UniformNumberGenerator(const double& raMin, const double& raMax, const double& decMin, const double& decMax);
  unsigned int& getNbObj(void);
  // DEFAULT CTOR
  HTMAsciiParser(void);

  // DEFAULT DTOR
  ~HTMAsciiParser(void);
 private:
  HTM* _htm;
  unsigned int nbObj;
};

}
