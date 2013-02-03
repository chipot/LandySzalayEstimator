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
#include "log.hh"

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

                this->_cache.emplace_back(nb1, nb2);
                //this->_htm->itemsToStore(nb1, nb2);
                //this->_htm->AddPoint(nb1, nb2);
                nbObj++;
            }
        }
        file.close();
    }
    else
        llog::fatal["HTMAsciiParser"]
            << "Can't Open File : " << filename << std::endl;
}

void HTMAsciiParser::Populate()
{
    for (auto const &p : this->_cache)
    {
        this->_htm->itemsToStore(p.first, p.second);
        this->_htm->AddPoint(p.first, p.second);
    }
}

unsigned int& HTMAsciiParser::getNbObj(void)
{
    return this->nbObj;
}

void    HTMAsciiParser::UniformNumberGenerator(const double& raMin, const double& raMax, const double& decMin, const double& decMax)
{
    unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 gen(seed1);
    std::uniform_real_distribution<double> unif1(raMin, raMax);
    std::uniform_real_distribution<double> unif2(decMin, decMax);
    auto random_ra = [&] () -> double {return unif1(gen);};
    auto random_dec = [&] () -> double {return unif2(gen);};

    for (unsigned int i = 0; i < nbObj; ++i)
    {
        const double ra = random_ra();
        const double dec = random_dec();
        this->_htm->AddPoint(ra, dec);
    }
}

HTMAsciiParser::HTMAsciiParser(HTM *htm)
    : _htm{htm}
    , nbObj(0)
{
}

HTMAsciiParser::~HTMAsciiParser()
{
}

}
