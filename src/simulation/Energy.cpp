#include "Energy.h"

char* Energy::model_t = "fem";
utils::Logger Energy::energy_logs("energy");
unsigned long Energy::ENERGY_EVALUATION_COUNTER(0);

Energy::Energy() {}

Energy::Energy(const Energy& orig) {}

Energy::~Energy() {}

double Energy::calcTotalEnergy(const std::vector<Shell>& shells, const Box& box)
{
    ENERGY_EVALUATION_COUNTER++;
    return 0.0;
}

double Energy::calcMembraneEnergy(const std::vector<Shell>& shells)
{
    return 0.0;
}

double Energy::calcStretchEnergy(const std::vector<Shell>& shells)
{
    return 0.0;
}

double Energy::calcBendingEnergy(const std::vector<Shell>& shells)
{
    return 0.0;
}

double Energy::calcOsmoticEnergy(const std::vector<Shell>& shells)
{
    double osmoticEnergy = 0.0;

    return osmoticEnergy;
}

double Energy::calcContactEnergy(const std::vector<Shell>& shells, const Box& box)
{
    return 0.0;
}

double Energy::calcShellBoxEnergy(const std::vector<Shell>& shells, const Box& box)
{
    return 0.0;
}

double Energy::calcShellShellEnergy(const std::vector<Shell>& shells, const Box& box)
{
    return 0.0;
}