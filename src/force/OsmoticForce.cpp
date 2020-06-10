#include "OsmoticForce.h"

void OsmoticForce::set_volume_flag(bool flag)
{
    volume_flag = flag;
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
    return volume_flag;
}

double OsmoticForce::epsilon = 0.0;
bool OsmoticForce::volume_flag = false;