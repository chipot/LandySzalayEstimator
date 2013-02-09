// C++ includes
#include <string>
#include <fstream>
#include <tbb/pipeline.h>
#include <tbb/concurrent_vector.h>
#include <boost/algorithm/string.hpp>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <cstdlib>

//BLINK Logservice
#include "log.hh"

//BLINK includes
#include "htmasciiparser.hpp"

namespace htm
{

tbb::concurrent_vector<std::pair<double, double>>
HTMAsciiParser::Parse(std::string& filename)
{
    std::ifstream   file(filename);
    //std::vector<std::pair<double, double>> points;
    tbb::concurrent_vector<std::pair<double, double>> points;

    struct stat ss;
    int fd = open(filename.c_str(), O_RDONLY);
    fstat(fd, &ss);
    void *map = mmap(nullptr, ss.st_size, PROT_WRITE | PROT_READ, MAP_PRIVATE, fd, 0);

    char *start = (char*)map;
    char *end = (char *)((intptr_t)map + ss.st_size);
    if (file)
    {
        tbb::parallel_pipeline(10000,

                               tbb::make_filter<void, char *>(
                               tbb::filter::serial,
                               [&] (tbb::flow_control &fc) -> char * {
                                 static char *it = start;

                                 char *line_begin = it;
                                 if (it == end)
                                 {
                                   fc.stop();
                                   return NULL;
                                 }
                                 while (*it != '\n') ++it;
                                 *it = '\0';
                                 ++it;
                                 return line_begin;
                               }) &

                               tbb::make_filter<char *, void>(
                               tbb::filter::parallel,
                               [&] (char *line) -> void {
                                 char *it = line;

                                 while (*it != ' ') ++it;
                                 *it = '\0';
                                 double ra = std::strtod(line, nullptr);
                                 ++it;
                                 double dec = std::strtod(it, nullptr);
                                 points.push_back(std::make_pair(ra, dec));
                                 return ;
                               })
        );
    }
    else
        llog::fatal["HTMAsciiParser"]
            << "Can't Open File : " << filename << std::endl;
    points.shrink_to_fit();
    return std::move(points);
}

}
