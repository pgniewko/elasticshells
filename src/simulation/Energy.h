#ifndef ENERGY_H
#define	ENERGY_H

#include "Environment.h"
#include "simulation/Box.h"
#include "Cell.h"
#include "simulation/DomainList.h"

class Energy
{
    public:
        Energy();
        Energy(const Energy& orig);
        virtual ~Energy();

        static double calcTotalEnergy(const std::vector<Cell>&, const Box&, const DomainList&);
        static double calcMembraneEnergy(const std::vector<Cell>&);
        static double calcOsmoticEnergy(const std::vector<Cell>&);
        static double calcContactEnergy(const std::vector<Cell>&, const Box&, const DomainList&);

        static double calcBendingEnergy(const std::vector<Cell>&);
        static double calcStretchEnergy(const std::vector<Cell>&);
        static double calcCellCellEnergy(const std::vector<Cell>&, const Box&, const DomainList&);
        static double calcCellBoxEnergy(const std::vector<Cell>&, const Box&);

        static void setModelType(char* _model_t)
        {
            model_t = _model_t;
            energy_logs << utils::LogLevel::INFO << "ENERGY MODEL HAS BEEN SET TO: " << model_t << "\n";
        }

        static ulong ENERGY_EVALUATION_COUNTER;

    private:

        static char* model_t;
        static utils::Logger energy_logs;
};

#endif	/* ENERGY_H */

