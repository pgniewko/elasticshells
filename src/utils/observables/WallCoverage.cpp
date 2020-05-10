#include "WallCoverage.h"

WallCoverage::WallCoverage(const char* name, const char* format) : Observer(name, format) {}

WallCoverage::WallCoverage(const WallCoverage& orig) : Observer(orig) {}

WallCoverage::~WallCoverage() {}

bool WallCoverage::is_touching_box(const Box& box, const Vector3D& vertex, 
        const double v_r, const double eps=0.0)
{
    double vx = vertex.x;
    double vy = vertex.y;
    double vz = vertex.z;

    double bx = box.get_x() - v_r - eps;
    double by = box.get_y() - v_r - eps;
    double bz = box.get_z() - v_r - eps;

    if (vx > bx || vx < -bx)
        return true;
    if (vy > by || vy < -by)
        return true;
    if (vz > bz || vz < -bz)
        return true;

    return false;   
}

bool WallCoverage::is_in_contact(const Box& box, const Shell& shell, const uint t_idx)
{
    if (box.pbc)
    {
        return false;
    }

    int idx1 = shell.triangles[t_idx].ia;
    int idx2 = shell.triangles[t_idx].ib;
    int idx3 = shell.triangles[t_idx].ic;

    bool a = is_touching_box(box, shell.vertices[idx1].r_c, shell.get_vertex_size(), 1e-5);
    bool b = is_touching_box(box, shell.vertices[idx2].r_c, shell.get_vertex_size(), 1e-5);
    bool c = is_touching_box(box, shell.vertices[idx3].r_c, shell.get_vertex_size(), 1e-5);

    if (a && b && c)
    {
        return true;
    }

    return false;
}


double WallCoverage::contact_area(const Box& box, const Shell& shell)
{ 
    uint t_nums = shell.get_number_triangles();
    double cell_wall_contact_area = 0.0;

    Vector3D cm = shell.get_cm();

    for (uint tid = 0; tid < t_nums; tid++)
    {
        bool iic = is_in_contact(box, shell, tid);
        if (iic)
        {
            cell_wall_contact_area += shell.triangles[tid].area(shell.vertices, cm, params[0]);
        }
    }

    return cell_wall_contact_area;
}


double WallCoverage::observe(const Box& box, const std::vector<Shell>& shells)
{
    if (box.pbc)
    {
        return 0.0;
    }
 
    double box_area = box.get_area(params[1]);

    uint cells_number = shells.size();
    double coverage = 0.0;
 
    for (uint i = 0; i < cells_number; i++)
    {
        coverage += contact_area(box, shells[i]);
    }

    coverage /= box_area;
    return coverage;
}

DerivedRegister<WallCoverage> WallCoverage::reg("WallCoverage");
