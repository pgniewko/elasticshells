#include "Tinker.h" 

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
            cell.vertices[cell.numberV] = Vertex(xtmp, ytmp, ztmp);
            cell.vertices[cell.numberV].setId(cell.numberV);
            cell.numberV++;
        }

        if ( Tinker::isUnique(vectors, i->b) )
        {
            vectors.push_back(i->b);
            xtmp = i->b.x;
            ytmp = i->b.y;
            ztmp = i->b.z;
            cell.vertices[cell.numberV] = Vertex(xtmp, ytmp, ztmp);
            cell.vertices[cell.numberV].setId(cell.numberV);
            cell.numberV++;
        }

        if ( Tinker::isUnique(vectors, i->c) )
        {
            vectors.push_back(i->c);
            xtmp = i->c.x;
            ytmp = i->c.y;
            ztmp = i->c.z;
            cell.vertices[cell.numberV] = Vertex(xtmp, ytmp, ztmp);
            cell.vertices[cell.numberV].setId(cell.numberV);
            cell.numberV++;
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
        cell.triangles[cell.numberT] = VertexTriangle(vrxt);
        cell.triangles[cell.numberT].setId(cell.numberT);
        cell.numberT++;
    }
}

int Tinker::getVertex(Cell& cell, const Vector3D& v)
{
    for (int i = 0; i < cell.numberV; i++)
    {
        if (cell.vertices[i].xyz.x == v.x && cell.vertices[i].xyz.y == v.y && cell.vertices[i].xyz.z == v.z)
        {
            return cell.vertices[i].getId();
        }
    }

    return -1;
}

void Tinker::constructTopology(Cell& cell)
{
    for (int i = 0; i < cell.numberT; i++)
    {
        int aid = cell.triangles[i].ia;
        int bid = cell.triangles[i].ib;
        int cid = cell.triangles[i].ic;
        Vector3D ab = cell.vertices[aid].xyz - cell.vertices[bid].xyz;
        Vector3D ac = cell.vertices[aid].xyz - cell.vertices[cid].xyz;
        Vector3D bc = cell.vertices[bid].xyz - cell.vertices[cid].xyz;
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

void Tinker::grow(Cell& cell)
{
    int vertexId = getRandomVertex(cell);
    
    
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
    else if(vertPos == 1)
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

    
    Vector3D newcoor = 0.5 * (cell.vertices[vert1].xyz + cell.vertices[vert2].xyz);
    cell.vertices[cell.numberV] = Vertex(newcoor.x, newcoor.y, newcoor.z);
    cell.vertices[cell.numberV].setId(cell.numberV);
    cell.vertices[cell.numberV].setMass(cell.vertices[vertexId].getMass());
    cell.vertices[cell.numberV].setVisc(cell.vertices[vertexId].getVisc());
    int newid = cell.numberV;
    cell.numberV++;
    
    
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
    
    Vector3D ab;// = vertices[aid].xyz - vertices[bid].xyz;
    ab = cell.vertices[newid].xyz - cell.vertices[vert1].xyz;
    //vertices[newid].addNeighbor(vert1, ab.length());
    //vertices[vert1].addNeighbor(newid, ab.length());
    
    ab = cell.vertices[newid].xyz - cell.vertices[vert2].xyz;
    //vertices[newid].addNeighbor(vert2, ab.length());
    //vertices[vert2].addNeighbor(newid, ab.length());
    
    ab = cell.vertices[newid].xyz - cell.vertices[vertexId].xyz;
    //vertices[newid].addNeighbor(vertexId, ab.length());
    //vertices[vertexId].addNeighbor(newid, ab.length());
    
    ab = cell.vertices[newid].xyz - cell.vertices[vert3].xyz;
    //vertices[newid].addNeighbor(vert3, ab.length() );
    //vertices[vert3].addNeighbor(newid, ab.length() );
    
    // =========
    ab = cell.vertices[newid].xyz - cell.vertices[vert1].xyz;
    cell.vertices[newid].addNeighbor(vert1, ab.length() );
    cell.vertices[vert1].addNeighbor(newid, ab.length() );
    
    ab = cell.vertices[newid].xyz - cell.vertices[vert2].xyz;
    cell.vertices[newid].addNeighbor(vert2, ab.length() );
    cell.vertices[vert2].addNeighbor(newid, ab.length() );
    
    ab = cell.vertices[newid].xyz - cell.vertices[vertexId].xyz;
    cell.vertices[newid].addNeighbor(vertexId, ab.length() );
    cell.vertices[vertexId].addNeighbor(newid, ab.length() );
    
    ab = cell.vertices[newid].xyz - cell.vertices[vert3].xyz;
    cell.vertices[newid].addNeighbor(vert3, ab.length() );
    cell.vertices[vert3].addNeighbor(newid, ab.length() );
    
    // ======
    //vertices[newid].addNeighbor(vert1, 0.5 * r0av);
    //vertices[vert1].addNeighbor(newid, 0.5 * r0av);
    
    //vertices[newid].addNeighbor(vert2, 0.5 * r0av);
    //vertices[vert2].addNeighbor(newid, 0.5 * r0av);
    
    //vertices[newid].addNeighbor(vertexId, 0.5 * r0av);
    //vertices[vertexId].addNeighbor(newid, 0.5 * r0av);
    
    //vertices[newid].addNeighbor(vert3, 0.5 * r0av);
    //vertices[vert3].addNeighbor(newid, 0.5 * r0av);
    // ==========
    

    
    // REMOVE TRIANGLES FROM THE VERTICES RECORD
    cell.vertices[vert2].removeTriangle(triangleId);
    cell.vertices[vert2].removeTriangle(secondTriangleId);
    
    cell.triangles[triangleId].subsVertex(vert2, newid); 
    cell.triangles[secondTriangleId].subsVertex(vert2, newid);
    cell.vertices[newid].addTriangle(triangleId);
    cell.vertices[newid].addTriangle(secondTriangleId);
    
    // NEW TRIANGLE 1
    VertexTriangle newTri1(vertexId, vert2, newid);
    cell.triangles[cell.numberT] = VertexTriangle(newTri1);
    cell.triangles[cell.numberT].setId(cell.numberT);
    int tid1 = cell.triangles[cell.numberT].myindex;
    cell.vertices[vertexId].addTriangle(tid1);
    cell.vertices[vert2].addTriangle(tid1);
    cell.vertices[newid].addTriangle(tid1);
    cell.numberT++;
    
    
    // NEW TRIANGLE 2
    VertexTriangle newTri2(newid, vert2, vert3);
    cell.triangles[cell.numberT] = VertexTriangle(newTri2);
    cell.triangles[cell.numberT].setId(cell.numberT);
    int tid2 = cell.triangles[cell.numberT].myindex;
    cell.vertices[newid].addTriangle(tid2);
    cell.vertices[vert2].addTriangle(tid2);
    cell.vertices[vert3].addTriangle(tid2);
    cell.numberT++;
    
    cell.vcounter++;
    ///std::cout << "vcounter=" << vcounter << std::endl;

    
}

// CREATE LIST DEP ON CELL CYCLE STAGE
int Tinker::getRandomVertex(Cell& cell)
{
    int vertexId = -1;
    int vertexCounter = 0;
    if (cell.my_phase == cell_phase_t::C_G1)
    {
        for (int i = 0; i < cell.numberV; i++)
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
        for (int i = 0; i < cell.numberV; i++)
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
        std::cout << "ERROR a.k.a. BUG" << std::endl;
        exit(1);
    }
    
    // RANDOM
    vertexId = uniform(0, vertexCounter);
    return vidx[vertexId];
    
    //double spatialNb[cell.numberV];
    //for (int i = 0; i < cell.numberV; i++)
    //{
    //    spatialNb[i] = 0.0;
    //}
    
    //double ptot = 0.0;
    //double cutoff = cell.r0av;
    //for (int i = 0; i< cell.numberV; i++)
    //{
    //    for (int j = 0; j < cell.numberV; j++)
    //    {
    //        double r = (cell.vertices[i].xyz - cell.vertices[j].xyz).length() ;
    //        if (r <= cutoff && j!=i) 
    //        {
    //            spatialNb[i] += 1.0 / r;
                //spatialNb[j] += 1.0;// / r;
                //ptot += 1.0;// / r;
                //ptot += 1.0;// / r;
    //        }
    //    }
    //}
    
    //for (int i = 0; i < cell.numberV; i++)
    //{
    //    spatialNb[i] = std::max(spatialNb[i], 1.0);
    //}
    
    //for (int i = 0; i < cell.numberV; i++)
    //{
    //    ptot += 1.0 / spatialNb[i];
    //    std::cout << "spatialNb[" << i << "]=" << spatialNb[i] << std::endl;
    //}
    //std::cout << "ptot=" << ptot << std::endl;
    
    
    //double sumcheck = 0.0;
    //for (int i = 0; i < cell.numberV; i++)
    //{
    //    sumcheck += (1. / spatialNb[i]) / (ptot);
    //}
    
    //double randn = uniform(0, 1.0);
    //double fracsum = 0.0;
   
    //for (int i = 0; i < cell.numberV; i++)
    //{
    //    fracsum += (1. / spatialNb[i]) / (ptot);
    //    if (randn < fracsum)
    //    {
    //        vertexId = i;
    //        break;
    //    }
    //}

    //return vertexId;
}

void Tinker::bud(Cell& cell)
{
    return;
}

void Tinker::divide(Cell& cell)
{
    return;
}

int Tinker::vidx[MAX_V];