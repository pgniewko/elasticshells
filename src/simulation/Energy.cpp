#include "Energy.h"

char* Energy::model_t = "fem";
utils::Logger Energy::energy_logs("energy");
unsigned long Energy::ENERGY_EVALUATION_COUNTER(0);

Energy::Energy() {}

Energy::Energy(const Energy& orig) {}

Energy::~Energy() {}

double Energy::calc_total_energy(const std::vector<Shell>& shells, const Box& box)
{
    ENERGY_EVALUATION_COUNTER++;
    return 0.0;
}

double Energy::calc_membrane_energy(const std::vector<Shell>& shells)
{
    return 0.0;
}

double Energy::calc_stretch_energy(const std::vector<Shell>& shells)
{
    return 0.0;
}

double Energy::calc_bending_energy(const std::vector<Shell>& shells)
{
    return 0.0;
}

double Energy::calc_osmotic_energy(const std::vector<Shell>& shells)
{
    double osmoticEnergy = 0.0;

    return osmoticEnergy;
}

double Energy::calc_contact_energy(const std::vector<Shell>& shells, const Box& box)
{
    return 0.0;
}

double Energy::calc_shell_box_energy(const std::vector<Shell>& shells, const Box& box)
{
    return 0.0;
}

double Energy::calc_shell_shell_energy(const std::vector<Shell>& shells, const Box& box)
{
    return 0.0;
}