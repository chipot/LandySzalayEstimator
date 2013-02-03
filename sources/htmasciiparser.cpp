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
#include "htmasciiparser.hpp"

namespace htm
{

void HTMAsciiParser::Parse(std::string& filename)
{
    std::ifstream   file(filename);
    std::string line(""), str1(""), str2("");
    size_t  pos1, pos2;
    nbObj = 0;

    if (file)
    {
        while (getline(file, line))
        {
            pos1 = line.find_first_of(" ", 0);
            str1 = line.substr(0, pos1);
            pos2 = line.find_first_of(" ", (pos1 + 1));
            str2 = line.substr(pos1+1, pos2-pos1-1);
            if ((str1 != " ") && (str1 != "\n") && (str1 != ""))
            {
                const double nb1 = strtod(str1.c_str(), NULL);
                const double nb2 = strtod(str2.c_str(), NULL);

                this->_htm->itemsToStore(nb1, nb2);
                this->_htm->AddPoint(nb1, nb2);
                nbObj++;
            }
        }
        file.close();
    }
    else
        LS_ADDMSG(LogService::FATAL, "HTMAsciiParser", "Can't Open File : " + filename);
}

unsigned int& HTMAsciiParser::getNbObj(void)
{
    return this->nbObj;
}

void    HTMAsciiParser::UniformNumberGenerator(const double& raMin, const double& raMax, const double& decMin, const double& decMax)
{
    std::random_device rd;
    unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 gen(seed1);
    std::uniform_real_distribution<double> unif1(raMin, raMax);
    std::uniform_real_distribution<double> unif2(decMin, decMax);
    for (unsigned int i = 0; i < nbObj; ++i)
    {
        const double ra = unif1(gen);
        const double dec = unif2(gen);
        this->_htm->AddPoint(ra, dec);
    }
}

HTMAsciiParser::HTMAsciiParser()
: nbObj(0)
{
    this->_htm = HTM::GetInstance();
}

HTMAsciiParser::~HTMAsciiParser()
{
}

}
