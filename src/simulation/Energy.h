#ifndef ENERGY_H
#define	ENERGY_H

#include "Environment.h"
#include "Shell.h"
#include "simulation/Box.h"

class Energy
{
    public:
        Energy();
        Energy(const Energy& orig);
        virtual ~Energy();

        static double calcTotalEnergy(const std::vector<Shell>&, const Box&);
        static double calcMembraneEnergy(const std::vector<Shell>&);
        static double calcOsmoticEnergy(const std::vector<Shell>&);
        static double calcContactEnergy(const std::vector<Shell>&, const Box&);

        static double calcBendingEnergy(const std::vector<Shell>&);
        static double calcStretchEnergy(const std::vector<Shell>&);
        static double calcShellShellEnergy(const std::vector<Shell>&, const Box&);
        static double calcShellBoxEnergy(const std::vector<Shell>&, const Box&);

        static unsigned long ENERGY_EVALUATION_COUNTER;

    private:

        static char* model_t;
        static utils::Logger energy_logs;
};

#endif	/* ENERGY_H */

