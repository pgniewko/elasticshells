#include "Pressure.h"

utils::Logger Pressure::sp_log("pressure");

Pressure::Pressure(const char* name, const char* format) : Observer(name, format) {}

Pressure::Pressure(const Pressure& orig) : Observer(orig) {}

Pressure::~Pressure() {}

void Pressure::set_params(const int num, std::vector<std::string> args_)
{
    d_param = strtod(args_[ num + 0 ].c_str(), NULL);
}

double Pressure::observe(const Box& box, const std::vector<Shell>& shells)
{    
    double pressure = 0.0;

    if (box.pbc && info_not_printed)
    {
        sp_log << utils::LogLevel::INFO << "PRESSURE FOR PBC IS NOT IMPLEMENTED AT THE MOMENT" << "\n";
        info_not_printed = false;
        return 0.0;
    }

    double totalForce = total_force(box, shells);
    double area = box.get_area(d_param);
    pressure = totalForce / area;

    return pressure;
}

double Pressure::total_force(const Box& box, const std::vector<Shell>& shells)
{
    if (box.pbc)
    {
        return 0.0;
    }
    if (image_not_created)
    {
        create_shells_image(box, shells);
    }
    copy_shells_data(box, shells);
    
    Vector3D wall_yz(0, 0, 0);
    Vector3D wall_xz(0, 0, 0);
    Vector3D wall_xy(0, 0, 0);
    Vector3D force_collector(0, 0, 0);
    Vector3D djk(0, 0, 0);

    double sgnx, sgny, sgnz;
    double bsx = box.get_x();
    double bsy = box.get_y();
    double bsz = box.get_z();

    double h;
    double e_eff;
    double rv = shells[0].get_vertex_size();
    double nu = shells[0].get_nu();
    double E = shells[0].get_E();
    double nub = box.get_nu();
    double Eb = box.get_E();

    double x, y, z;
    double tot_force = 0.0;

    uint n = xyz.size() / 3;
    Vector3D vertex;

    for (uint i = 0; i < n; i++)
    {
        force_collector *= 0;
        x = xyz[3 * i + 0];
        y = xyz[3 * i + 1];
        z = xyz[3 * i + 2];
        vertex = Vector3D(x, y, z);

        //////////////
        // WALL XY  //
        //////////////
        sgnx = SIGN(vertex.x);
        wall_yz.x = sgnx * bsx;
        wall_yz.y = vertex.y;
        wall_yz.z = vertex.z;
        djk = vertex - wall_yz;


        h = rv - djk.length();
        e_eff = (1 - nu * nu) / E + (1 - nub * nub) / Eb;
        e_eff = 1.0 / e_eff;

        if (h > 0)
        {
            tot_force += constants::d4_3 * e_eff * pow(rv, 0.5) * pow(h, 1.5);
            //force_collector += fmagn * (djk / djk.length());
        }

        //////////////
        // WALL XZ  //
        //////////////
        sgny = SIGN(vertex.y);
        wall_xz.x = vertex.x;
        wall_xz.y = sgny * bsy;
        wall_xz.z = vertex.z;
        djk = vertex - wall_xz;

        h = rv - djk.length();
        e_eff = (1 - nu * nu) / E + (1 - nub * nub) / Eb;
        e_eff = 1.0 / e_eff;

        if (h > 0)
        {
            tot_force += constants::d4_3 * e_eff * pow(rv, 0.5) * pow(h, 1.5);
            //force_collector += fmagn * (djk / djk.length());
        }

        //////////////
        // WALL XY  //
        //////////////
        sgnz = SIGN(vertex.z);
        wall_xy.x = vertex.x;
        wall_xy.y = vertex.y;
        wall_xy.z = sgnz * bsz;
        djk = vertex - wall_xy;

        h = rv - djk.length();
        e_eff = (1 - nu * nu) / E + (1 - nub * nub) / Eb;
        e_eff = 1.0 / e_eff;

        if (h > 0)
        {
            tot_force += constants::d4_3 * e_eff * pow(rv, 0.5) * pow(h, 1.5);
            //force_collector += fmagn * (djk / djk.length());
        }

        //tot_force += force_collector.length();

        //forces[3 * i + 0] += force_collector.x;
        //forces[3 * i + 1] += force_collector.y;
        //forces[3 * i + 2] += force_collector.z;
    }

    return tot_force;
}

DerivedRegister<Pressure> Pressure::reg("Pressure");