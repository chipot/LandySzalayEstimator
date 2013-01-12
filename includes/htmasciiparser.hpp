#pragma once

// C++ includes
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <random>
#include <chrono>

// C includes
#include <math.h>
#include <time.h>

//BLINK Logservice
#include "logservice.hpp"

//BLINK includes
#include "htm.hpp"

using namespace ICoDF;
using namespace ICoDF_HTM;

namespace ICoDF_HTM
{
class HTM;
class HTMAsciiParser //: public IHTMParser
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
