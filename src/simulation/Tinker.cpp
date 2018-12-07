#include "Tinker.h"

utils::Logger Tinker::tinker_log("tinker");

Tinker::Tinker() {}

Tinker::Tinker(const Tinker& orig) {}

Tinker::~Tinker() {}

void Tinker::constructVertices(Shell& shell, std::list<Triangle>& tris)
{
    std::list<Vector3D> vectors;
    double xtmp, ytmp, ztmp;

    for (std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i)
    {
        if ( Tinker::isUnique(vectors, i->a) )
        {
            vectors.push_back(i->a);
            xtmp = i->a.x;
            ytmp = i->a.y;
            ztmp = i->a.z;
            shell.vertices.push_back(Vertex(xtmp, ytmp, ztmp));
            shell.vertices[shell.number_v].set_id(shell.number_v);
            shell.number_v++;
        }

        if ( Tinker::isUnique(vectors, i->b) )
        {
            vectors.push_back(i->b);
            xtmp = i->b.x;
            ytmp = i->b.y;
            ztmp = i->b.z;
            shell.vertices.push_back(Vertex(xtmp, ytmp, ztmp));
            shell.vertices[shell.number_v].set_id(shell.number_v);
            shell.number_v++;
        }

        if ( Tinker::isUnique(vectors, i->c) )
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

void Tinker::constructVTriangles(Shell& shell, std::list<Triangle>& tris)
{
    for (std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i)
    {
        int va = getVertex(shell, i->a);
        int vb = getVertex(shell, i->b);
        int vc = getVertex(shell, i->c);
        Element vrxt(va, vb, vc);
        shell.triangles.push_back(Element(vrxt));
        shell.triangles[shell.number_t].setId(shell.number_t);
        shell.number_t++;
    }
}

int Tinker::getVertex(Shell& cell, const Vector3D& v, double e)
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

void Tinker::constructTopology(Shell& shell)
{
    for (int i = 0; i < shell.number_t; i++)
    {
        int aid = shell.triangles[i].ia;
        int bid = shell.triangles[i].ib;
        int cid = shell.triangles[i].ic;
        int tid = shell.triangles[i].my_id;

        shell.vertices[aid].addNeighbor(bid);
        shell.vertices[aid].addNeighbor(cid);
        shell.vertices[bid].addNeighbor(aid);
        shell.vertices[bid].addNeighbor(cid);
        shell.vertices[cid].addNeighbor(aid);
        shell.vertices[cid].addNeighbor(bid);
        shell.vertices[aid].addTriangle(tid);
        shell.vertices[bid].addTriangle(tid);
        shell.vertices[cid].addTriangle(tid);
    }
}

void Tinker::constructBSprings(Shell& shell)
{

    if (shell.getNumberVertices() == 1)
    {
        shell.number_h = 0;
        return;
    }

    for (int x3 = 0; x3 < shell.number_v; x3++)
    {
        for (int j = 0; j < shell.vertices[x3].vertex_degree; j++)
        {
            std::vector<int> common_verts;
            int x4 = shell.vertices[x3].bondedVerts[j];

            for (int k = 0; k < shell.vertices[x3].vertex_degree; k++)
            {
                for (int l = 0; l < shell.vertices[x4].vertex_degree; l++)
                {
                    if (x4 != shell.vertices[x3].bondedVerts[k] && x3 != shell.vertices[x4].bondedVerts[l])
                    {
                        if (shell.vertices[x3].bondedVerts[k] == shell.vertices[x4].bondedVerts[l])
                        {
                            common_verts.push_back(shell.vertices[x3].bondedVerts[k]);
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

                    if ( isBSpringUnique(x1_, x2_, x3_, x4_, shell) )
                    {
                        //cell.bhinges[cell.number_s] = BendingHinge(x1_, x2_, x3_, x4_);
                        shell.hinges.push_back(Hinge(x1_, x2_, x3_, x4_));
                        shell.hinges[shell.number_h].setId(shell.number_h);
                        shell.number_h++;
                    }
                }
            }
        }
    }
}

bool Tinker::isUnique(std::list<Vector3D>& vlist, Vector3D& v, double e)
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

bool Tinker::isBSpringUnique(int x1, int x2, int x3, int x4, Shell& cell)
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