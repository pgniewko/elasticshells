#include "Tinker.h"

utils::Logger Tinker::tinker_log("tinker");

Tinker::Tinker() {}

Tinker::Tinker(const Tinker& orig) {}

Tinker::~Tinker() {}

void Tinker::construct_vertices(Shell& shell, std::list<Triangle>& tris)
{
    std::list<Vector3D> vectors;
    double xtmp, ytmp, ztmp;

    for (std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i)
    {
        if ( Tinker::is_unique(vectors, i->a) )
        {
            vectors.push_back(i->a);
            xtmp = i->a.x;
            ytmp = i->a.y;
            ztmp = i->a.z;
            shell.vertices.push_back(Vertex(xtmp, ytmp, ztmp));
            shell.vertices[shell.number_v].set_id(shell.number_v);
            shell.number_v++;
        }

        if ( Tinker::is_unique(vectors, i->b) )
        {
            vectors.push_back(i->b);
            xtmp = i->b.x;
            ytmp = i->b.y;
            ztmp = i->b.z;
            shell.vertices.push_back(Vertex(xtmp, ytmp, ztmp));
            shell.vertices[shell.number_v].set_id(shell.number_v);
            shell.number_v++;
        }

        if ( Tinker::is_unique(vectors, i->c) )
        {
            vectors.push_back(i->c);
            xtmp = i->c.x;
            ytmp = i->c.y;
            ztmp = i->c.z;
            shell.vertices.push_back(Vertex(xtmp, ytmp, ztmp));
            shell.vertices[shell.number_v].set_id(shell.number_v);
            shell.number_v++;
        }
    }
}

void Tinker::construct_elements(Shell& shell, std::list<Triangle>& tris)
{
    for (std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i)
    {
        int va = get_vertex(shell, i->a);
        int vb = get_vertex(shell, i->b);
        int vc = get_vertex(shell, i->c);
        Element vrxt(va, vb, vc);
        shell.triangles.push_back(Element(vrxt));
        shell.triangles[shell.number_t].set_id(shell.number_t);
        shell.number_t++;
    }
}

int Tinker::get_vertex(Shell& cell, const Vector3D& v, double e)
{
    for (int i = 0; i < cell.number_v; i++)
    {
        if ( fabs(cell.vertices[i].r_c.x - v.x) < e && fabs(cell.vertices[i].r_c.y - v.y) < e && fabs(cell.vertices[i].r_c.z - v.z) < e )
        {
            return cell.vertices[i].get_id();
        }
    }
    return -1;
}

void Tinker::construct_topology(Shell& shell)
{
    for (int i = 0; i < shell.number_t; i++)
    {
        int aid = shell.triangles[i].ia;
        int bid = shell.triangles[i].ib;
        int cid = shell.triangles[i].ic;
        int tid = shell.triangles[i].my_id;

        shell.vertices[aid].add_neighbor(bid);
        shell.vertices[aid].add_neighbor(cid);
        shell.vertices[bid].add_neighbor(aid);
        shell.vertices[bid].add_neighbor(cid);
        shell.vertices[cid].add_neighbor(aid);
        shell.vertices[cid].add_neighbor(bid);
        shell.vertices[aid].add_element(tid);
        shell.vertices[bid].add_element(tid);
        shell.vertices[cid].add_element(tid);
    }
}

void Tinker::construct_hinges(Shell& shell)
{

    if (shell.get_number_vertices() == 1)
    {
        shell.number_h = 0;
        return;
    }

    for (int x3 = 0; x3 < shell.number_v; x3++)
    {
        for (int j = 0; j < shell.vertices[x3].vertex_degree; j++)
        {
            std::vector<int> common_verts;
            int x4 = shell.vertices[x3].bonded_vertices[j];

            for (int k = 0; k < shell.vertices[x3].vertex_degree; k++)
            {
                for (int l = 0; l < shell.vertices[x4].vertex_degree; l++)
                {
                    if (x4 != shell.vertices[x3].bonded_vertices[k] && x3 != shell.vertices[x4].bonded_vertices[l])
                    {
                        if (shell.vertices[x3].bonded_vertices[k] == shell.vertices[x4].bonded_vertices[l])
                        {
                            common_verts.push_back(shell.vertices[x3].bonded_vertices[k]);
                        }
                    }
                }
            }

            int x3_ = std::min(x3, x4);
            int x4_ = std::max(x3, x4);

            for (uint ix = 0; ix < common_verts.size(); ix++)
            {
                for (uint iy = ix + 1; iy < common_verts.size(); iy++)
                {
                    int x1_ = std::min(common_verts[ix], common_verts[iy]);
                    int x2_ = std::max(common_verts[ix], common_verts[iy]);

                    if ( is_hinge_unique(x1_, x2_, x3_, x4_, shell) )
                    {
                        shell.hinges.push_back(Hinge(x1_, x2_, x3_, x4_));
                        shell.hinges[shell.number_h].setId(shell.number_h);
                        shell.number_h++;
                    }
                }
            }
        }
    }
}

bool Tinker::is_unique(std::list<Vector3D>& vlist, Vector3D& v, double e)
{
    for (std::list<Vector3D>::iterator i = vlist.begin(); i != vlist.end(); ++i)
    {
        if ( fabs(i->x - v.x) < e && fabs(i->y - v.y) < e && fabs(i->z - v.z) < e )
        {
            return false;
        }
    }

    return true;
}

bool Tinker::is_hinge_unique(int x1, int x2, int x3, int x4, Shell& cell)
{
    Hinge bs_tmp(x1, x2, x3, x4);

    for (int i = 0; i < cell.number_h; i++)
    {
        if (bs_tmp == cell.hinges[i])
        {
            return false;
        }
    }

    return true;
}