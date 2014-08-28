#include "HookeanForce.h"

HookeanForce::HookeanForce() {}

HookeanForce::HookeanForce(const HookeanForce& orig) {}

HookeanForce::~HookeanForce() {}

 static Vector3D HookeanForce::calcForce(const Vector3D& va, const Vector3D& vb, const double R0, const double gamma)
    {
        Vector3D dR = vb - va;
        double R = dR.length();
        Vector3D f = gamma * dR * (1 - R0 / R);
        return f;
    }