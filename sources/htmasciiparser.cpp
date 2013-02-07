// C++ includes
#include <fstream>

//BLINK Logservice
#include "log.hh"

//BLINK includes
#include "htmasciiparser.hpp"

namespace htm
{

std::vector<std::pair<double, double>>
HTMAsciiParser::Parse(std::string& filename)
{
    std::ifstream   file(filename);
    std::vector<std::pair<double, double>> points;

    if (file)
    {
        while (file.eof() == false)
        {
            double nb1;
            double nb2;

            file >> nb1 >> nb2;

            points.emplace_back(nb1, nb2);
        }
        file.close();
    }
    else
        llog::fatal["HTMAsciiParser"]
            << "Can't Open File : " << filename << std::endl;
    return std::move(points);
}

HTMAsciiParser::HTMAsciiParser(HTM *htm)
    : _htm{htm}
{
}

HTMAsciiParser::~HTMAsciiParser()
{
}

}
