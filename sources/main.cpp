#include <iostream>
#include <limits>
#include <unistd.h>
#include <sstream>

#include "htm.hpp"
#include "htmasciiparser.hpp"
#include "log.hh"

int main(int ac, char **av)
{
    llog::notice["main"] << "BLINK::HTM test main" << std::endl;
    if (ac == 4)
    {
        int loop = 100;
        std::string file(av[1]);
        htm::HTM *htm = htm::HTM::GetInstance();
        auto* parser = new htm::HTMAsciiParser;

        std::istringstream stm, stm2;
        stm.str(av[2]);
        double radius = 0;
        stm >> radius;
        stm2.str(av[3]);
        double delta = 0;
        stm2 >> delta;

        llog::notice["main"] << "Computing Normal Catalog..." << std::endl;;
        htm->CreateOctahedron();
        parser->Parse(file);
        htm->CreateHTM();
        llog::notice["main"] << "HTM Created for Normal Catalog" << std::endl;;
        unsigned int nn = htm->TwoPointsCorrelation(radius, delta);
        std::stringstream tmp;
        llog::notice["main"] << "Two Point Correlation have been computed for the Normal Catalog [" << nn << "] pairs" << std::endl; 

        double raMin = htm->getMinRa();
        double raMax = htm->getMaxRa();
        double decMin = htm->getMinDec();
        double decMax = htm->getMaxDec();
        llog::notice["main"] << "Computing Mean values for Random and Hybrid Catalog on " << loop << " loops using " << parser->getNbObj() << " random objects..." << std::endl;

        unsigned int rr = 0;
        unsigned int nr = 0;
        for (int i = 0; i != loop; ++i)
        {
            llog::notice["main"] << "Computing loop " << i + 1 << " on " << loop << std::endl;
            htm->DeleteOctahedron();
            htm->CreateOctahedron();
            parser->UniformNumberGenerator(raMin, raMax, decMin, decMax);
            htm->CreateHTM();

            unsigned int currentRR = htm->TwoPointsCorrelation(radius, delta);
            rr += currentRR;
            llog::notice["main"] << "Two Point Correlation have been computed for the Random Catalog [" << currentRR << "] mean [" << (rr / (i + 1)) << "]" << std::endl;

            parser->Parse(file);
            htm->CreateHTM();
            unsigned int currentNR = htm->TwoPointsCorrelation(radius, delta);
            nr += currentNR;
            llog::notice["main"] << "Two Point Correlation have been computed for the Hybrid Catalog [" << currentNR << "] mean [" << (nr / (i + 1)) << "]" << std::endl;
        }

        llog::notice["main"] << "...Done !" << std::endl;

        rr /= loop;
        nr /= loop;
        double estimator = 0;
        estimator = 2 * nr;
        estimator = nn - estimator;
        estimator = estimator - rr;
        estimator /= rr;
        llog::notice["main"] << "Landy Szalay Estimator for current catalog returns : " << estimator << std::endl;

        htm->Delete();
        delete parser;
    }
    else
        llog::fatal["main"] << "Usage : ./LandySzalayEstimator <catalog_file> <radius> <delta>" << std::endl;

    llog::notice["main"] << "Exiting..." << std::endl;
    return 0;
}
