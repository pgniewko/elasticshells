#include "Tinker.h"

utils::Logger Tinker::tinker_log("tinker");

int Tinker::vidx[MAX_V];

Tinker::Tinker() {}

Tinker::Tinker(const Tinker& orig) {}

Tinker::~Tinker() {}

void Tinker::constructVertices(Cell& cell, std::list<Triangle>& tris)
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
    }
}

bool Tinker::isUnique(std::list<Vector3D>& vlist, Vector3D& v)
{
    for (std::list<Vector3D>::iterator i = vlist.begin(); i != vlist.end(); ++i)
    {
        if (i->x == v.x && i->y == v.y && i->z == v.z)
        {
            return false;
        }
    }

    return true;
}


void Tinker::constructVTriangles(Cell& cell, std::list<Triangle>& tris)
{
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

int Tinker::getVertex(Cell& cell, const Vector3D& v)
{
    for (int i = 0; i < cell.number_v; i++)
    {
        if (cell.vertices[i].r_c.x == v.x && cell.vertices[i].r_c.y == v.y && cell.vertices[i].r_c.z == v.z)
        {
            return cell.vertices[i].getId();
        }
    }

    return -1;
}

void Tinker::constructTopology(Cell& cell)
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
        int tid = cell.triangles[i].myindex;
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

void Tinker::constructBSprings(Cell& cell)
{

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
                    //std::cout << cell.vertices[x3].bondedVerts[k] << " " << cell.vertices[x4].bondedVerts[l] << std::endl;
                    if (x4 != cell.vertices[x3].bondedVerts[k] && x3 != cell.vertices[x4].bondedVerts[l])
                    {
                        if (cell.vertices[x3].bondedVerts[k] == cell.vertices[x4].bondedVerts[l])
                        {
                            //std::cout << "ADD=" << cell.vertices[x3].bondedVerts[k] << std::endl;
                            common_verts.push_back(cell.vertices[x3].bondedVerts[k]);
                        }
                    }
                }
            }    
        
        
            int x3_ = std::min(x3, x4);
            int x4_ = std::max(x3, x4);
        
            //std::cout << "x3="<<x3 << " x4=" << x4 << " common_verts.size()="<<common_verts.size() << std::endl;
            for (uint ix = 0; ix < common_verts.size(); ix++)
            {
                for (uint iy = ix+1; iy < common_verts.size(); iy++)
                {
                    int x1_ = std::min(common_verts[ix], common_verts[iy]);
                    int x2_ = std::max(common_verts[ix], common_verts[iy]);
                    //std::cout << "x1_="<<x1_ << " x2_=" << x2_ << std::endl;
                
                    if( isBSpringUnique(x1_, x2_, x3_, x4_, cell) )
                    {
                        cell.bhinges[cell.number_s] = BendingSpring(x1_, x2_, x3_, x4_);
                        cell.number_s++;
                    }
                }
            }
        }
        //std::cout << "************" << std::endl;
        
        
    }
}

bool Tinker::isBSpringUnique(int x1, int x2, int x3, int x4, Cell& cell)
{
    BendingSpring bs_tmp(x1, x2, x3, x4);
    
    for (int i = 0; i < cell.number_s; i++)
    {
        if (bs_tmp == cell.bhinges[i])
        {
            return false;
        }
    }
    
    return true;
}

void Tinker::grow(Cell& cell)
{
    //int vertexId = getRandomVertex(cell);
    //int vertexId = getLonelyVertex(cell);
    int vertexId = getOldestVertex(cell);
    int triangle_num = uniform(0, cell.vertices[vertexId].numTris);
    int triangleId = cell.vertices[vertexId].bondedTris[triangle_num];
    int vertPos = -1;
    int vert1 = -1;
    int vert2 = -1;

    if (cell.triangles[triangleId].ia == vertexId)
    {
        vertPos = 0;
    }
    else if (cell.triangles[triangleId].ib == vertexId)
    {
        vertPos = 1;
    }
    else
    {
        vertPos = 2;
    }

    if (vertPos == 0)
    {
        vert1 = cell.triangles[triangleId].ib;
        vert2 = cell.triangles[triangleId].ic;
    }
    else if (vertPos == 1)
    {
        vert1 = cell.triangles[triangleId].ia;
        vert2 = cell.triangles[triangleId].ic;
    }
    else
    {
        vert1 = cell.triangles[triangleId].ia;
        vert2 = cell.triangles[triangleId].ib;
    }

    int secondTriangleId = -1;

    for (int i = 0; i < cell.vertices[vert1].numTris; i++)
    {
        for (int j = 0; j < cell.vertices[vert2].numTris; j++)
        {
            int t1 = cell.vertices[vert1].bondedTris[i];
            int t2 = cell.vertices[vert2].bondedTris[j];

            if (t1 == t2 && t1 != triangleId)
            {
                secondTriangleId = t1;
            }
        }
    }

    Vector3D newcoor = 0.5 * (cell.vertices[vert1].r_c + cell.vertices[vert2].r_c);
    cell.vertices[cell.number_v] = Vertex(newcoor.x, newcoor.y, newcoor.z);
    cell.vertices[cell.number_v].setId(cell.number_v);
    cell.vertices[cell.number_v].setVisc(cell.vertices[vertexId].getVisc());
    int newid = cell.number_v;
    cell.number_v++;
    int vert3 = -1;

    if (cell.triangles[secondTriangleId].ia == vert1)
    {
        if (cell.triangles[secondTriangleId].ib == vert2)
        {
            vert3 = cell.triangles[secondTriangleId].ic;
        }
        else
        {
            vert3 = cell.triangles[secondTriangleId].ib;
        }
    }

    if (cell.triangles[secondTriangleId].ib == vert1)
    {
        if (cell.triangles[secondTriangleId].ia == vert2)
        {
            vert3 = cell.triangles[secondTriangleId].ic;
        }
        else
        {
            vert3 = cell.triangles[secondTriangleId].ia;
        }
    }

    if (cell.triangles[secondTriangleId].ic == vert1)
    {
        if (cell.triangles[secondTriangleId].ia == vert2)
        {
            vert3 = cell.triangles[secondTriangleId].ib;
        }
        else
        {
            vert3 = cell.triangles[secondTriangleId].ia;
        }
    }

    //REMOVE OLD BONDS
    cell.vertices[vert1].removeNeighbor(vert2);
    cell.vertices[vert2].removeNeighbor(vert1);
    // CREATE NEW BONDS
    Vector3D ab;
    // /*
    ab = cell.vertices[newid].r_c - cell.vertices[vert1].r_c;
    cell.vertices[newid].addNeighbor(vert1, ab.length() );
    cell.vertices[vert1].addNeighbor(newid, ab.length() );
    ab = cell.vertices[newid].r_c - cell.vertices[vert2].r_c;
    cell.vertices[newid].addNeighbor(vert2, ab.length() );
    cell.vertices[vert2].addNeighbor(newid, ab.length() );
    ab = cell.vertices[newid].r_c - cell.vertices[vertexId].r_c;
    cell.vertices[newid].addNeighbor(vertexId, ab.length() );
    cell.vertices[vertexId].addNeighbor(newid, ab.length() );
    ab = cell.vertices[newid].r_c - cell.vertices[vert3].r_c;
    cell.vertices[newid].addNeighbor(vert3, ab.length() );
    cell.vertices[vert3].addNeighbor(newid, ab.length() );
    //*/
    /*
    cell.vertices[newid].addNeighbor(vert1, cell.r0av+uniform(-0.05,0.05));
    cell.vertices[vert1].addNeighbor(newid, cell.r0av+uniform(-0.05,0.05));


    cell.vertices[newid].addNeighbor(vert2, cell.r0av+uniform(-0.05,0.05));
    cell.vertices[vert2].addNeighbor(newid, cell.r0av+uniform(-0.05,0.05));


    cell.vertices[newid].addNeighbor(vertexId, cell.r0av+uniform(-0.05,0.05));
    cell.vertices[vertexId].addNeighbor(newid, cell.r0av+uniform(-0.05,0.05));


    cell.vertices[newid].addNeighbor(vert3, cell.r0av+uniform(-0.05,0.05));
    cell.vertices[vert3].addNeighbor(newid, cell.r0av+uniform(-0.05,0.05));
    */
    // REMOVE TRIANGLES FROM THE VERTICES RECORD
    cell.vertices[vert2].removeTriangle(triangleId);
    cell.vertices[vert2].removeTriangle(secondTriangleId);
    cell.triangles[triangleId].subsVertex(vert2, newid);
    cell.triangles[secondTriangleId].subsVertex(vert2, newid);
    cell.vertices[newid].addTriangle(triangleId);
    cell.vertices[newid].addTriangle(secondTriangleId);
    // NEW TRIANGLE 1
    VertexTriangle newTri1(vertexId, vert2, newid);
    cell.triangles[cell.number_t] = VertexTriangle(newTri1);
    cell.triangles[cell.number_t].setId(cell.number_t);
    int tid1 = cell.triangles[cell.number_t].myindex;
    cell.vertices[vertexId].addTriangle(tid1);
    cell.vertices[vert2].addTriangle(tid1);
    cell.vertices[newid].addTriangle(tid1);
    cell.number_t++;
    // NEW TRIANGLE 2
    VertexTriangle newTri2(newid, vert2, vert3);
    cell.triangles[cell.number_t] = VertexTriangle(newTri2);
    cell.triangles[cell.number_t].setId(cell.number_t);
    int tid2 = cell.triangles[cell.number_t].myindex;
    cell.vertices[newid].addTriangle(tid2);
    cell.vertices[vert2].addTriangle(tid2);
    cell.vertices[vert3].addTriangle(tid2);
    cell.number_t++;
    //cell.vcounter++;
    ///std::cout << "vcounter=" << vcounter << std::endl;
}

// CREATE LIST DEP ON CELL CYCLE STAGE
int Tinker::getRandomVertex(Cell& cell)
{
    int vertexCounter = validVertList(cell);
    int vertexId = uniform(0, vertexCounter);
    return vidx[vertexId];
}

int Tinker::getLonelyVertex(Cell& cell)
{
    int vertexCounter = validVertList(cell);
    double spatialNb[MAX_V];

    for (int i = 0; i < vertexCounter; i++)
    {
        spatialNb[i] = 0.0;
    }

    double ptot = 0.0;
    double cutoff = 2.0 * cell.params.vertex_r;
    int I, J;

    for (int i = 0; i < vertexCounter; i++)
    {
        for (int j = 0; j < vertexCounter; j++)
        {
            I = vidx[i];
            J = vidx[j];
            double r = (cell.vertices[J].r_c - cell.vertices[I].r_c).length() ;

            if (r <= cutoff && J != I)
            {
                spatialNb[i] += 1.0 / r;
            }
        }
    }

    for (int i = 0; i < vertexCounter; i++)
    {
        spatialNb[i] = std::max(spatialNb[i], 1.0);
    }

    for (int i = 0; i < vertexCounter; i++)
    {
        ptot += 1.0 / spatialNb[i];
    }

    double sumcheck = 0.0;

    for (int i = 0; i < cell.number_v; i++)
    {
        sumcheck += (1. / spatialNb[i]) / (ptot);
    }

    double randn = uniform(0, 1.0);
    double fracsum = 0.0;
    int vertexId = -1;

    for (int i = 0; i < vertexCounter; i++)
    {
        fracsum += (1. / spatialNb[i]) / (ptot);

        if (randn < fracsum)
        {
            vertexId = i;
            break;
        }
    }

    return vertexId;
}

int Tinker::getOldestVertex(Cell& cell)
{
    int vertexCounter = validVertList(cell);
    double vertAge[MAX_V];

    for (int i = 0; i < vertexCounter; i++)
    {
        vertAge[i] = 0.0;
    }

    double ptot = 0.0;
    int I;
//    int J;

    for (int i = 0; i < vertexCounter; i++)
    {
        I = vidx[i];
        vertAge[i] += cell.vertices[I].gtimer;
    }

    for (int i = 0; i < vertexCounter; i++)
    {
        ptot += vertAge[i];
    }

    if (ptot == 0.0)
    {
        std::cerr << "ptot=0.0. Bug - FIX IT !" << std::endl;
        exit(EXIT_FAILURE);
    }

    double randn = uniform(0, 1.0);
    double fracsum = 0.0;
    int vertexId = -1;

    for (int i = 0; i < vertexCounter; i++)
    {
        fracsum += vertAge[i] / ptot;

        if (randn < fracsum)
        {
            vertexId = i;
            break;
        }
    }

    //std::cout << vertAge[vertexId] << std::endl;
    cell.vertices[vertexId].voidTime();
    return vertexId;
}

int Tinker::validVertList(Cell& cell)
{
    int vertexCounter = 0;

    if (cell.my_phase == cell_phase_t::C_G1)
    {
        for (int i = 0; i < cell.number_v; i++)
        {
            if (cell.vertices[i].my_type == vertex_t::MOTHER)
            {
                vidx[vertexCounter] = i;
                vertexCounter++;
            }
        }
    }
    else if (cell.my_phase == cell_phase_t::C_SG2)
    {
        for (int i = 0; i < cell.number_v; i++)
        {
            if (cell.vertices[i].my_type == vertex_t::BUD)
            {
                vidx[vertexCounter] = i;
                vertexCounter++;
            }
        }
    }
    else
    {
        std::cout << "ERROR. Cell not is G1 or G2 phases." << std::endl;
        exit(EXIT_FAILURE);
    }

    return vertexCounter;
}

void Tinker::bud(Cell& cell)
{
    return;
}

void Tinker::divide(Cell& cell)
{
    return;
}