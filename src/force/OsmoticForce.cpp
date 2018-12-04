#include "OsmoticForce.h"

void OsmoticForce::setVolumeFlag(bool flag)
{
    volumeFlag = flag;
}

void OsmoticForce::setEpsilon(double eps)
{
    epsilon = eps;
}

double OsmoticForce::getEpsilon()
{
    return epsilon;
}

bool OsmoticForce::getFlag()
{
    return volumeFlag;
}

double OsmoticForce::epsilon = 0.0;
bool OsmoticForce::volumeFlag = false;