#include <iostream>
#include <limits>
#include <unistd.h>
#include <sstream>

#include "htm.hpp"
#include "htmasciiparser.hpp"
#include "log.hh"

static unsigned int loop = 100;

int main(int ac, char **av)
{
    llog::notice["main"] << "BLINK::HTM test main" << std::endl;
    if (ac == 4)
    {
        std::string file(av[1]);
        htm::HTM *htm = htm::HTM::GetInstance();
        htm::HTMAsciiParser parser(htm);

        double radius = std::stod(av[2]);
        double delta = std::stod(av[3]);

        llog::notice["main"] << "Computing Normal Catalog..." << std::endl;;
        htm->CreateOctahedron();
        parser.Parse(file);
        parser.PopulateHTM();
        htm->CreateHTM();

        llog::notice["main"] << "HTM Created for Normal Catalog" << std::endl;;
        unsigned int nn = htm->TwoPointsCorrelation(radius, delta);
        llog::notice["main"] << "Two Point Correlation have been computed for the Normal Catalog [" << nn << "] pairs" << std::endl; 

        double raMin = htm->getMinRa();
        double raMax = htm->getMaxRa();
        double decMin = htm->getMinDec();
        double decMax = htm->getMaxDec();
        llog::notice["main"] << "Computing Mean values for Random and Hybrid Catalog on " << loop << " loops using " << parser.getNbObj() << " random objects..." << std::endl;

        unsigned int rr = 0;
        unsigned int nr = 0;
        for (unsigned int i = 0; i != loop; ++i)
        {
            llog::notice["main"] << "Computing loop " << i + 1 << " on " << loop << std::endl;
            htm->DeleteOctahedron();
            htm->CreateOctahedron();
            parser.UniformNumberGenerator(raMin, raMax, decMin, decMax);
            htm->CreateHTM();

            unsigned int currentRR = htm->TwoPointsCorrelation(radius, delta);
            rr += currentRR;
            llog::notice["main"] << "Two Point Correlation have been computed for the Random Catalog [" << currentRR << "] mean [" << (rr / (i + 1)) << "]" << std::endl;

            parser.PopulateHTM();
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
    }
    else
        llog::fatal["main"] << "Usage : ./LandySzalayEstimator <catalog_file> <radius> <delta>" << std::endl;

    llog::notice["main"] << "Exiting..." << std::endl;
    return 0;
}
