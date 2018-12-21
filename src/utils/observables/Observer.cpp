#include "Observer.h"

Observer::Observer(const char* name, const char* format) :
    observer_name(name), output_format(format), i_param(0), d_param(0.0) {}

Observer::Observer(const Observer& orig) :
    observer_name(orig.observer_name), output_format(orig.output_format), i_param(orig.i_param), d_param(orig.d_param) {};

Observer::~Observer() {}

const char* Observer::getFormat()
{
    return output_format.c_str();
}

const char* Observer::getName()
{
    return observer_name.c_str();
}

void Observer::create_shells_image(const Box& box, const std::vector<Shell>& shells)
{
    int vertex_counter = 0;

    for (uint i = 0; i < shells.size(); i++)
    {
        turgors.push_back( 0.0 );
        double x_, y_, z_;

        for (int j = 0; j < shells[i].get_number_vertices(); j++)
        {
            x_ = shells[i].vertices[j].r_c.x;
            y_ = shells[i].vertices[j].r_c.y;
            z_ = shells[i].vertices[j].r_c.z;

            xyz.push_back(x_);
            xyz.push_back(y_);
            xyz.push_back(z_);

            forces.push_back(0.0);
            forces.push_back(0.0);
            forces.push_back(0.0);

            object_map vm(i, j);
            vs_map.push_back(vm);

            inv_vs_map[vm] = vertex_counter;

            std::vector<int> bonds;
            graph.push_back(bonds);

            vertex_counter++;
        }
    }

    for (uint i = 0; i < shells.size(); i++)
    {
        for (int j = 0; j < shells[i].get_number_vertices(); j++)
        {
            object_map vm_ij(i, j);
            int ij_id = inv_vs_map[vm_ij];

            for (int k = 0; k < shells[i].get_number_vertices(); k++)
            {
                if (k != j && shells[i].vertices[j].is_neighbor(k))
                {
                    object_map vm_ik(i, k);
                    int ik_id = inv_vs_map[vm_ik];
                    graph[ij_id].push_back(ik_id);
                }
            }

        }
    }


    int element_counter = 0;

    for (uint i = 0; i < shells.size(); i++)
    {
        int ia, ib, ic;
        int ia_mapped, ib_mapped, ic_mapped;

        for (int j = 0; j < shells[i].get_number_triangles(); j++)
        {
            element el_;
            ia = shells[i].triangles[j].ia;
            ib = shells[i].triangles[j].ib;
            ic = shells[i].triangles[j].ic;

            object_map vm_a(i, ia);
            object_map vm_b(i, ib);
            object_map vm_c(i, ic);

            ia_mapped = inv_vs_map[vm_a];
            ib_mapped = inv_vs_map[vm_b];
            ic_mapped = inv_vs_map[vm_c];

            el_.ia = ia_mapped;
            el_.ib = ib_mapped;
            el_.ic = ic_mapped;

            el_.an[0] = shells[i].triangles[j].an[0];
            el_.an[1] = shells[i].triangles[j].an[1];
            el_.an[2] = shells[i].triangles[j].an[2];

            el_.L2[0] = shells[i].triangles[j].L2[0];
            el_.L2[1] = shells[i].triangles[j].L2[1];
            el_.L2[2] = shells[i].triangles[j].L2[2];

            el_.ki[0] = shells[i].triangles[j].ki[0];
            el_.ki[1] = shells[i].triangles[j].ki[1];
            el_.ki[2] = shells[i].triangles[j].ki[2];

            el_.ci[0] = shells[i].triangles[j].ci[0];
            el_.ci[1] = shells[i].triangles[j].ci[1];
            el_.ci[2] = shells[i].triangles[j].ci[2];

            elements.push_back(el_);

            object_map vm(i, j);
            ts_map.push_back(vm);
            inv_ts_map[vm] = element_counter;

            element_counter++;
        }
    }

    int hinge_counter = 0;

    for (uint i = 0; i < shells.size(); i++)
    {
        int x1, x2, x3, x4;
        int x1_mapped, x2_mapped, x3_mapped, x4_mapped;

        for (int j = 0; j < shells[i].get_number_hinges(); j++)
        {

            hinge h_;
            x1 = shells[i].hinges[j].x1;
            x2 = shells[i].hinges[j].x2;
            x3 = shells[i].hinges[j].x3;
            x4 = shells[i].hinges[j].x4;

            object_map vm_a(i, x1);
            object_map vm_b(i, x2);
            object_map vm_c(i, x3);
            object_map vm_d(i, x4);

            x1_mapped = inv_vs_map[vm_a];
            x2_mapped = inv_vs_map[vm_b];
            x3_mapped = inv_vs_map[vm_c];
            x4_mapped = inv_vs_map[vm_d];

            h_.v1 = x1_mapped;
            h_.v2 = x2_mapped;
            h_.v3 = x3_mapped;
            h_.v4 = x4_mapped;

            h_.D = shells[i].hinges[j].D;
            h_.theta0 = shells[i].hinges[j].theta0;
            h_.sinTheta0 = shells[i].hinges[j].sinTheta0;

            hinges.push_back(h_);


            object_map vm(i, j);
            hs_map.push_back(vm);
            inv_hs_map[vm] = hinge_counter;
            hinge_counter++;
        }
    }
    image_not_created = false;
}

void Observer::copy_shells_data(const Box& box, const std::vector<Shell>& shells)
{
    int vertex_no = 0;

    for (uint i = 0; i < shells.size(); i++)
    {
        double x_, y_, z_;

        for (int j = 0; j < shells[i].get_number_vertices(); j++)
        {
            x_ = shells[i].vertices[j].r_c.x;
            y_ = shells[i].vertices[j].r_c.y;
            z_ = shells[i].vertices[j].r_c.z;

            object_map vm(i, j);

            vertex_no =  inv_vs_map[ vm ];

            xyz[3 * vertex_no + 0] = x_;
            xyz[3 * vertex_no + 1] = y_;
            xyz[3 * vertex_no + 2] = z_;
        }
    }

    for (uint i = 0; i < shells.size(); i++)
    {
        turgors[i] = shells[i].get_turgor();
    }

    fc.set_dl_dims(-box.get_x(), box.get_x(), 0);
    fc.set_dl_dims(-box.get_y(), box.get_y(), 1);
    fc.set_dl_dims(-box.get_z(), box.get_z(), 2);

    contacts = fc.contacts_list(xyz, graph, vs_map, (int) shells.size(), shells[0].get_vertex_size());
}

bool Observer::is_in_contact(int i, int j)
{
    std::vector<int> line = contacts[i];

    for (uint k = 0; k < line.size(); k++)
        if (line[k] == j)
        {
            return true;
        }

    return false;
}

ObserverFactory::map_type* ObserverFactory::map = NULL;