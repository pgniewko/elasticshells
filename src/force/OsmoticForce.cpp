#include "OsmoticForce.h"

void OsmoticForce::set_volume_flag(bool flag)
{
    volumeFlag = flag;
}

void OsmoticForce::set_epsilon(double eps)
{
    epsilon = eps;
}

double OsmoticForce::get_epsilon()
{
    return epsilon;
}

bool OsmoticForce::get_flag()
{
    return volumeFlag;
}

double OsmoticForce::epsilon = 0.0;
bool OsmoticForce::volumeFlag = false;