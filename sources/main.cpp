#include <iostream>
#include <limits>
#include <unistd.h>
#include <sstream>

#include "logservice.hpp"
#include "htm.hpp"
#include "htmasciiparser.hpp"

int main(int ac, char **av)
{
    // SETTING UP LOG
    LogService::GetInstance()->SetConfiguration(LogService::LS_PRINT_ON_COUT);
    LS_ADDMSG(LogService::NOTICE, "main", "BLINK::HTM test main");
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

        LS_ADDMSG(LogService::NOTICE, "main", "Computing Normal Catalog...");
        htm->CreateOctahedron();
        parser->Parse(file);
        htm->CreateHTM();
        LS_ADDMSG(LogService::NOTICE, "main", "HTM Created for Normal Catalog");
        unsigned int nn = htm->TwoPointsCorrelation(radius, delta);
        std::stringstream tmp;
        tmp.str("");
        tmp << "Two Point Correlation have been computed for the Normal Catalog [" << nn << "] pairs"; 
        LS_ADDMSG(LogService::NOTICE, "main", tmp.str());

        double raMin = htm->getMinRa();
        double raMax = htm->getMaxRa();
        double decMin = htm->getMinDec();
        double decMax = htm->getMaxDec();
        tmp.str("");
        tmp << "Computing Mean values for Random and Hybrid Catalog on " << loop << " loops using " << parser->getNbObj() << " random objects...";
        LS_ADDMSG(LogService::NOTICE, "main", tmp.str());

        unsigned int rr = 0;
        unsigned int nr = 0;
        for (int i = 0; i != loop; ++i)
        {
            tmp.str("");
            tmp << "Computing loop " << i + 1 << " on " << loop; 
            LS_ADDMSG(LogService::NOTICE, "main", tmp.str());
            htm->DeleteOctahedron();
            htm->CreateOctahedron();
            parser->UniformNumberGenerator(raMin, raMax, decMin, decMax);
            htm->CreateHTM();

            unsigned int currentRR = htm->TwoPointsCorrelation(radius, delta);
            rr += currentRR;
            tmp.str("");
            tmp << "Two POint Correlation have been computed for the Random Catalog [" << currentRR << "] mean [" << (rr / (i + 1)) << "]";

            LS_ADDMSG(LogService::NOTICE, "main", tmp.str());

            parser->Parse(file);
            htm->CreateHTM();
            unsigned int currentNR = htm->TwoPointsCorrelation(radius, delta);
            nr += currentNR;
            tmp.str("");
            tmp << "Two POint Correlation have been computed for the Hybrid Catalog [" << currentNR << "] mean [" << (nr / (i + 1)) << "]";
            LS_ADDMSG(LogService::NOTICE, "main", tmp.str());
        }

        LS_ADDMSG(LogService::NOTICE, "main", "...Done !");

        rr /= loop;
        nr /= loop;
        double estimator = 0;
        estimator = 2 * nr;
        estimator = nn - estimator;
        estimator = estimator - rr;
        estimator /= rr;
        tmp.str("");
        tmp << "Landy Szalay Estimator for current catalog returns : " << estimator;
        LS_ADDMSG(LogService::NOTICE, "main", tmp.str());

        htm->Delete();
        delete parser;
    }
    else
        LS_ADDMSG(LogService::FATAL, "main", "Usage : ./LandySzalayEstimator <catalog_file> <radius> <delta>");

    LS_ADDMSG(LogService::NOTICE, "main", "Exiting...");
    LogService::Delete();
    return 0;
}
