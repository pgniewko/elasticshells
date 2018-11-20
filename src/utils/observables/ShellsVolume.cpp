#include "ShellsVolume.h"

ShellsVolume::ShellsVolume (const char* name, const char* format) : Observer(name, format) {}

ShellsVolume::ShellsVolume(const ShellsVolume& orig) : Observer(orig) {}

ShellsVolume::~ShellsVolume() {}

void ShellsVolume::set_params(const int num, std::vector<std::string> args_)
{
    d_param = strtod(args_[ num + 0 ].c_str(), NULL);
}

double ShellsVolume::observe(const Box& box, std::vector<Shell>& shells, const DomainList& dl)
{
    return calcShellsVolume(shells);
}

double ShellsVolume::calcShellsVolume(std::vector<Shell>& shells)
{
    double shellsvolume = 0.0;

    for (uint i = 0; i < shells.size(); i++)
    {
        shellsvolume += shells[i].calcVolume(d_param);
    }

    if ( shells[0].getNumberVertices() == 1 )
    {
        double overlapvolume = 0.0;

        for (uint i = 0; i < shells.size(); i++)
        {
            for (uint j = i + 1; j < shells.size(); j++)
            {
                overlapvolume += sphereSphereIntersection(shells[i], shells[j]);
            }
        }

        shellsvolume -= overlapvolume;
    }

    return shellsvolume;
}

double ShellsVolume::sphereSphereIntersection(const Shell& c1, const Shell& c2)
{

    double r1 = c1.getInitR();
    double r2 = c2.getInitR();
    Vector3D cm1 = c1.getCm();
    Vector3D cm2 = c2.getCm();

    double d = (cm1 - cm2).length();

    double V = 0.0;

    double rmin = std::min(r1, r2);

    if ( d > (r1 + r2) )
    {
        V = 0.0;
    }
    else if ( d <= std::abs(r1 - r2) )
    {
        V = 4.0 * constants::pi * rmin * rmin * rmin / 3.0;
    }
    else
    {
        V  = (r1 + r2 - d) * (r1 + r2 - d) * constants::pi;
        V *= (d * d + 2 * d * (r1 + r2) - 3 * (r1 * r1 + r2 * r2) + 6 * r1 * r2);
        V /= (12 * d);
    }

    return V;
}

DerivedRegister<ShellsVolume> ShellsVolume::reg("ShellsVolume");