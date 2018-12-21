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

        static double calc_total_energy(const std::vector<Shell>&, const Box&);
        static double calc_membrane_energy(const std::vector<Shell>&);
        static double calc_osmotic_energy(const std::vector<Shell>&);
        static double calc_contact_energy(const std::vector<Shell>&, const Box&);

        static double calc_bending_energy(const std::vector<Shell>&);
        static double calc_stretch_energy(const std::vector<Shell>&);
        static double calc_shell_shell_energy(const std::vector<Shell>&, const Box&);
        static double calc_shell_box_energy(const std::vector<Shell>&, const Box&);

        static unsigned long ENERGY_EVALUATION_COUNTER;

    private:

        static char* model_t;
        static utils::Logger energy_logs;
};

#endif	/* ENERGY_H */

