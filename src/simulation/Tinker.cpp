#include "Tinker.h"

utils::Logger Tinker::tinker_log("tinker");

Tinker::Tinker() {}

Tinker::Tinker(const Tinker& orig) {}

Tinker::~Tinker() {}

void Tinker::constructVertices(Shell& cell, std::list<Triangle>& tris)
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
            cell.vertices[cell.number_v] = Vertex(xtmp, ytmp, ztmp);
            cell.vertices[cell.number_v].setId(cell.number_v);
            cell.number_v++;
        }

        if ( Tinker::isUnique(vectors, i->b) )
        {
            vectors.push_back(i->b);
            xtmp = i->b.x;
            ytmp = i->b.y;
            ztmp = i->b.z;
            cell.vertices[cell.number_v] = Vertex(xtmp, ytmp, ztmp);
            cell.vertices[cell.number_v].setId(cell.number_v);
            cell.number_v++;
        }

        if ( Tinker::isUnique(vectors, i->c) )
        {
            vectors.push_back(i->c);
            xtmp = i->c.x;
            ytmp = i->c.y;
            ztmp = i->c.z;
            cell.vertices[cell.number_v] = Vertex(xtmp, ytmp, ztmp);
            cell.vertices[cell.number_v].setId(cell.number_v);
            cell.number_v++;
        }

        if (cell.number_v > MAX_V)
        {
            tinker_log << utils::LogLevel::CRITICAL  << "The number of vertices is larger than allowed in Environment.h\n";
            tinker_log << utils::LogLevel::CRITICAL  << "The simulation will be terminated ! \n";
            exit(EXIT_SUCCESS);
        }
    }
}

void Tinker::constructVTriangles(Shell& cell, std::list<Triangle>& tris)
{

    if (tris.size() > MAX_T)
    {
        tinker_log << utils::LogLevel::CRITICAL  << "The number of triangles is larger than allowed in Environment.h\n";
        tinker_log << utils::LogLevel::CRITICAL  << "The simulation will be terminated ! \n";
        exit(EXIT_SUCCESS);
    }

    for (std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i)
    {
        int va = getVertex(cell, i->a);
        int vb = getVertex(cell, i->b);
        int vc = getVertex(cell, i->c);
        VertexTriangle vrxt(va, vb, vc);
        cell.triangles[cell.number_t] = VertexTriangle(vrxt);
        cell.triangles[cell.number_t].setId(cell.number_t);
        cell.number_t++;
    }
}

int Tinker::getVertex(Shell& cell, const Vector3D& v, double e)
{
    for (int i = 0; i < cell.number_v; i++)
    {
        if ( fabs(cell.vertices[i].r_c.x - v.x) < e && fabs(cell.vertices[i].r_c.y - v.y) < e && fabs(cell.vertices[i].r_c.z - v.z) < e )
        {
            return cell.vertices[i].getId();
        }
    }

    //std::cout << " v = " << v<< std::endl;
    return -1;
}

void Tinker::constructTopology(Shell& cell)
{
    for (int i = 0; i < cell.number_t; i++)
    {
        int aid = cell.triangles[i].ia;
        int bid = cell.triangles[i].ib;
        int cid = cell.triangles[i].ic;
        Vector3D ab = cell.vertices[aid].r_c - cell.vertices[bid].r_c;
        Vector3D ac = cell.vertices[aid].r_c - cell.vertices[cid].r_c;
        Vector3D bc = cell.vertices[bid].r_c - cell.vertices[cid].r_c;
        double abl = ab.length();
        double acl = ac.length();
        double bcl = bc.length();
        int tid = cell.triangles[i].myid;
        cell.vertices[aid].addNeighbor(bid, abl);
        cell.vertices[aid].addNeighbor(cid, acl);
        cell.vertices[bid].addNeighbor(aid, abl);
        cell.vertices[bid].addNeighbor(cid, bcl);
        cell.vertices[cid].addNeighbor(aid, acl);
        cell.vertices[cid].addNeighbor(bid, bcl);
        cell.vertices[aid].addTriangle(tid);
        cell.vertices[bid].addTriangle(tid);
        cell.vertices[cid].addTriangle(tid);
    }
}

void Tinker::constructBSprings(Shell& cell)
{

    if (cell.getNumberVertices() == 1)
    {
        cell.number_s = 0;
        return;
    }

    for (int x3 = 0; x3 < cell.number_v; x3++)
    {
        for (int j = 0; j < cell.vertices[x3].numBonded; j++)
        {
            std::vector<int> common_verts;
            int x4 = cell.vertices[x3].bondedVerts[j];

            for (int k = 0; k < cell.vertices[x3].numBonded; k++)
            {
                for (int l = 0; l < cell.vertices[x4].numBonded; l++)
                {
                    if (x4 != cell.vertices[x3].bondedVerts[k] && x3 != cell.vertices[x4].bondedVerts[l])
                    {
                        if (cell.vertices[x3].bondedVerts[k] == cell.vertices[x4].bondedVerts[l])
                        {
                            common_verts.push_back(cell.vertices[x3].bondedVerts[k]);
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

                    if ( isBSpringUnique(x1_, x2_, x3_, x4_, cell) )
                    {
                        cell.bhinges[cell.number_s] = BendingHinge(x1_, x2_, x3_, x4_);
                        cell.bhinges[cell.number_s].setId(cell.number_s);
                        cell.number_s++;
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
    BendingHinge bs_tmp(x1, x2, x3, x4);

    for (int i = 0; i < cell.number_s; i++)
    {
        if (bs_tmp == cell.bhinges[i])
        {
            return false;
        }
    }

    return true;
}