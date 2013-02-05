#include <iostream>
#include <limits>
#include <unistd.h>
#include <sstream>
#include <iomanip>

#include "htm.hpp"
#include "htmasciiparser.hpp"
#include "log.hh"

static unsigned int loop = 100;

int main(int ac, char **av)
{
    llog::debug["main"] << "BLINK::HTM test main" << std::endl;
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

        llog::debug["main"] << "HTM Created for Normal Catalog" << std::endl;;
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
            llog::debug["main"] << "Computing loop " << i + 1 << " on " << loop << std::endl;
            htm->DeleteOctahedron(); // Oh my God
            htm->CreateOctahedron(); // So that's why...
            // This is stupid: This function add points into the HTM !!!!
            parser.UniformNumberGenerator(raMin, raMax, decMin, decMax);
            htm->CreateHTM(); // This one just pop points from a list filled by the previous call

            unsigned int currentRR = htm->TwoPointsCorrelation(radius, delta);
            rr += currentRR;
            llog::debug["main"] << "Two Point Correlation have been computed for the Random Catalog [" << currentRR << "] mean [" << (rr / (i + 1)) << "]" << std::endl;

            parser.PopulateHTM();
            htm->CreateHTM();
            unsigned int currentNR = htm->TwoPointsCorrelation(radius, delta);
            nr += currentNR;
            llog::debug["main"] << "Two Point Correlation have been computed for the Hybrid Catalog [" << currentNR << "] mean [" << (nr / (i + 1)) << "]" << std::endl;
            std::cout << "\r" << i;
        }
        std::cout << std::endl;

        llog::debug["main"] << "\n...Done !" << std::endl;

        // rr and nr are modified in the loop.. nn is defined by the Two P. C.
        // at the beginning
        rr /= loop;
        nr /= loop;

        double estimator = 0;
        {
            double nn_d = nn;
            double nr_d = nr;
            double rr_d = rr;

            estimator = (nn_d - (2 * nr_d) + rr_d) / rr_d;
        }

        llog::notice["main"] << "Landy Szalay Estimator for current catalog returns : "<< std::scientific << std::setprecision(15) << estimator << std::endl;

        htm->Delete();
    }
    else
        llog::fatal["main"] << "Usage : ./LandySzalayEstimator <catalog_file> <radius> <delta>" << std::endl;

    llog::debug["main"] << "Exiting..." << std::endl;
    return 0;
}
