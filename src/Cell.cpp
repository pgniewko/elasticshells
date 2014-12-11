#include "Cell.h"

Cell::Cell(int depth) :  cellId(-1), numberV(0), numberT(0), nRT(0), r0av(0), 
        vcounter(0), budVno(0), my_phase(cell_phase_t::C_G0)
{
    SimpleTriangulation sm(depth);
    std::list<Triangle> tris = sm.triangulate();
    constructVertices(tris);
    constructVTriangles(tris);
    constructTopology();
    calcR0av();
}

Cell::Cell(std::list<Triangle> tris) : cellId(-1), numberV(0), numberT(0), nRT(0), 
        r0av(0), vcounter(0), budVno(0), my_phase(cell_phase_t::C_G0)
{
    constructVertices(tris);
    constructVTriangles(tris);
    constructTopology();
    calcR0av();
}

Cell::Cell(const Cell& orig) : cm(orig.cm), vertices(orig.vertices), triangles(orig.triangles),
    cellId(orig.cellId), params(orig.params), numberV(orig.numberV), numberT(orig.numberT), nRT(orig.nRT), r0av(orig.r0av),
        vcounter(orig.vcounter), budVno(orig.budVno), my_phase(orig.my_phase)
{}

Cell::~Cell() {}

void Cell::constructVertices(std::list<Triangle> tris)
{
    std::list<Vector3D> vectors;
    double xtmp, ytmp, ztmp;

    for (std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i)
    {
        if ( isUnique(vectors, i->a) )
        {
            vectors.push_back(i->a);
            xtmp = i->a.x;
            ytmp = i->a.y;
            ztmp = i->a.z;
            vertices[numberV] = Vertex(xtmp, ytmp, ztmp);
            vertices[numberV].setId(numberV);
            numberV++;
        }

        if ( isUnique(vectors, i->b) )
        {
            vectors.push_back(i->b);
            xtmp = i->b.x;
            ytmp = i->b.y;
            ztmp = i->b.z;
            vertices[numberV] = Vertex(xtmp, ytmp, ztmp);
            vertices[numberV].setId(numberV);
            numberV++;
        }

        if ( isUnique(vectors, i->c) )
        {
            vectors.push_back(i->c);
            xtmp = i->c.x;
            ytmp = i->c.y;
            ztmp = i->c.z;
            vertices[numberV] = Vertex(xtmp, ytmp, ztmp);
            vertices[numberV].setId(numberV);
            numberV++;
        }
    }
}

bool Cell::isUnique(std::list<Vector3D>& vlist, const Vector3D& v)
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

void Cell::constructVTriangles(std::list<Triangle> tris)
{
    for (std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i)
    {
        int va = getVertex(i->a);
        int vb = getVertex(i->b);
        int vc = getVertex(i->c);
        VertexTriangle vrxt(va, vb, vc);
        triangles[numberT] = VertexTriangle(vrxt);
        triangles[numberT].setId(numberT);
        numberT++;
    }
}

int Cell::getVertex(const Vector3D v)
{
    for (int i = 0; i < numberV; i++)
    {
        if (vertices[i].xyz.x == v.x && vertices[i].xyz.y == v.y && vertices[i].xyz.z == v.z)
        {
            return vertices[i].getId();
        }
    }

    return -1;
}


void Cell::constructTopology()
{
    for (int i = 0; i < numberT; i++)
    {
        int aid = triangles[i].ia;
        int bid = triangles[i].ib;
        int cid = triangles[i].ic;
        Vector3D ab = vertices[aid].xyz - vertices[bid].xyz;
        Vector3D ac = vertices[aid].xyz - vertices[cid].xyz;
        Vector3D bc = vertices[bid].xyz - vertices[cid].xyz;
        double abl = ab.length();
        double acl = ac.length();
        double bcl = bc.length();
        int tid = triangles[i].myindex;
        vertices[aid].addNeighbor(bid, abl);
        vertices[aid].addNeighbor(cid, acl);
        vertices[bid].addNeighbor(aid, abl);
        vertices[bid].addNeighbor(cid, bcl);
        vertices[cid].addNeighbor(aid, acl);
        vertices[cid].addNeighbor(bid, bcl);
        vertices[aid].addTriangle(tid);
        vertices[bid].addTriangle(tid);
        vertices[cid].addTriangle(tid);
    }
}

void Cell::voidVerletLsit()
{
    for (int i = 0; i < numberV; i++)
    {
        vertices[i].numNbNeighs = 0;
    }
}

void Cell::builtVerletList(const Cell& other_cell, Box& box)
{
    Vector3D distance_jk;
    double r_cut = 2 * params.r_vertex;

    if (this->cellId != other_cell.cellId)
    {
        for (int j = 0; j < numberV; j++)
        {
            for (int k = 0; k < other_cell.numberV; k++ )
            {
                getDistance(distance_jk,  other_cell.vertices[k].xyz, vertices[j].xyz, box);

                if (distance_jk.length() <= r_cut * params.verletR)
                {
                    vertices[j].addNbNeighbor(k, other_cell.cellId);
                }
            }
        }
    }
    else
    {
        for (int j = 0; j < numberV; j++)
        {
            for (int k = 0; k < numberV; k++)
            {
                getDistance(distance_jk, vertices[k].xyz, vertices[j].xyz, box);

                if (j != k && !vertices[j].isNeighbor(k) && distance_jk.length() <= r_cut * params.verletR)
                {
                    vertices[j].addNbNeighbor(k, other_cell.cellId);
                }
            }
        }
    }
}

void Cell::builtNbList(std::vector<Cell>& cells, DomainList& domains, Box& box)
{
    int domainIdx;
    int domainn;
    Vector3D distance_ik;
    int vertIdx, cellIdx;

    double r_cut = 2 * params.r_vertex + EPSILON;
    
    for (int i = 0; i < numberV; i++)
    {
        domainIdx = vertices[i].domainIdx;
        domainIdx = domains.getDomainIndex(vertices[i]);

        for (int j = 0; j < domains.getNumberOfNeigh(domainIdx); j++)
        {
            domainn = domains.getDomainNeighbor(domainIdx, j);

            for (int k = 0; k < domains.getNumOfParticles(domainn); k++)
            {
                vertIdx = domains.getVertexIdx(domainn, k);
                cellIdx = domains.getCellIdx(domainn, k);

                if (this->cellId != cellIdx)
                {
                    getDistance(distance_ik, cells[cellIdx].vertices[vertIdx].xyz, vertices[i].xyz, box);

                    if (distance_ik.length() <= r_cut)
                    {
                        vertices[i].addNbNeighbor(vertIdx, cellIdx);
                    }
                }
                else
                {
                    getDistance(distance_ik, vertices[vertIdx].xyz, vertices[i].xyz, box);

                    if (i != vertIdx && !vertices[i].isNeighbor(vertIdx) && distance_ik.length() <= r_cut)
                    {
                        vertices[i].addNbNeighbor(vertIdx, cellIdx);
                    }
                }
            }
        }
    }
#ifdef TESTS
    for (int i = 0; i < numberV; i++)
    {
        vertices[i].sortNbList();
    }
#endif    
    
}

void Cell::calcBondedForces()
{
    calcHarmonicForces();
    calcOsmoticForces();
}

void Cell::calcHarmonicForces()
{
    double R0ij;
    int idxj;

    for (int i = 0; i < numberV; i++)
    {
        for (int j = 0; j < vertices[i].numBonded; j++)
        {
            R0ij = vertices[i].r0[j];
            idxj = vertices[i].bondedVerts[j];
            vertices[i].force += HookeanForce::calcForce(vertices[i].xyz, vertices[idxj].xyz, R0ij, params.gamma);
        }
    }
}

void Cell::calcOsmoticForces()
{
    calcCM();
    double cellVolume = calcVolume();
    int iva, ivb, ivc;
    double dp = params.dp;

    for (int i = 0; i < numberT; i++)
    {
        iva = triangles[i].ia;
        ivb = triangles[i].ib;
        ivc = triangles[i].ic;
        Vector3D fa = OsmoticForce::calcForce(vertices[iva].xyz, vertices[ivb].xyz, vertices[ivc].xyz, cm, nRT, cellVolume, dp);
        Vector3D fb = OsmoticForce::calcForce(vertices[ivb].xyz, vertices[ivc].xyz, vertices[iva].xyz, cm, nRT, cellVolume, dp);
        Vector3D fc = OsmoticForce::calcForce(vertices[ivc].xyz, vertices[iva].xyz, vertices[ivb].xyz, cm, nRT, cellVolume, dp);
        vertices[iva].force += fa;
        vertices[ivb].force += fb;
        vertices[ivc].force += fc;
    }
}

void Cell::calcNbForcesON2(const Cell& other_cell, Box& box)
{
    int ocellid = other_cell.cellId;
    Vector3D dij;

    for (int i = 0; i < numberV; i++)
    {
        for (int j = 0; j < other_cell.numberV; j++)
        {
            if (cellId != ocellid)
            {
                getDistance(dij, other_cell.vertices[j].xyz, vertices[i].xyz, box);
                vertices[i].force += HertzianRepulsion::calcForce(dij, params.r_vertex, params.r_vertex, params.ecc);
            }
            else
            {
                if (i != j && !vertices[i].isNeighbor(j))
                {
                    getDistance(dij, other_cell.vertices[j].xyz, vertices[i].xyz, box);
                    vertices[i].force += HertzianRepulsion::calcForce(dij, params.r_vertex, params.r_vertex, params.ecc);
                }
            }
        }
    }
}

void Cell::calcNbForcesVL(const Cell& other_cell, Box& box)
{
    int ocellid = other_cell.cellId;
    int vertid;
    Vector3D divix;

    for (int i = 0; i < numberV; i++)
    {
        for (int j = 0; j < vertices[i].numNbNeighs; j++)
        {
            if (vertices[i].nbCellsIdx[j] == ocellid)
            {
                vertid = vertices[i].nbVerts[j];
                getDistance(divix, other_cell.vertices[vertid].xyz, vertices[i].xyz, box);
                vertices[i].force += HertzianRepulsion::calcForce(divix, params.r_vertex, params.r_vertex, params.ecc);
            }
        }
    }
}

void Cell::calcBoxForces(Box& box)
{
    Vector3D wallYZ, wallXZ, wallXY;
    Vector3D dij;
    double sgnx, sgny, sgnz;
    double bsx = box.getX();
    double bsy = box.getY();
    double bsz = box.getZ();
    double ecw = box.ecw;

    for (int i = 0; i < numberV; i++)
    {
        sgnx = SIGN(vertices[i].xyz.x);
        wallYZ.x = sgnx * bsx;
        wallYZ.y = vertices[i].xyz.y;
        wallYZ.z = vertices[i].xyz.z;
        dij = vertices[i].xyz - wallYZ;
        vertices[i].force += HertzianRepulsion::calcForce(dij, params.r_vertex, ecw);
        
        sgny = SIGN(vertices[i].xyz.y);
        wallXZ.x = vertices[i].xyz.x;
        wallXZ.y = sgny * bsy;
        wallXZ.z = vertices[i].xyz.z;
        dij = vertices[i].xyz - wallXZ;
        vertices[i].force += HertzianRepulsion::calcForce(dij, params.r_vertex, ecw);
        
        sgnz = SIGN(vertices[i].xyz.z);
        wallXY.x = vertices[i].xyz.x;
        wallXY.y = vertices[i].xyz.y;
        wallXY.z = sgnz * bsz;
        dij = vertices[i].xyz - wallXY;
        vertices[i].force += HertzianRepulsion::calcForce(dij, params.r_vertex, ecw);
    }
}

void Cell::voidForces()
{
    for (int i = 0; i < numberV; i++)
    {
        vertices[i].voidForce();
    }
}

double Cell::calcSurfaceArea()
{
    double surface = 0.0;

    for (int i = 0; i < numberT; i++)
    {
        surface += triangles[i].area(vertices);
    }

    return surface;
}

double Cell::calcVolume()
{
    double volume = 0.0;
    int va, vb, vc;
    calcCM();

    for (int i = 0; i < numberT; i++)
    {
        va = triangles[i].ia;
        vb = triangles[i].ib;
        vc = triangles[i].ic;
        Tetrahedron tetr(vertices[va].xyz, vertices[vb].xyz, vertices[vc].xyz, cm);
        volume += tetr.volume();
    }

    return volume;
}

void Cell::calcCM()
{
    Vector3D tmp(0.0, 0.0, 0.0);
    double M = 0.0;
    double m;

    for (int i = 0; i < numberV; i++)
    {
        m = vertices[i].getMass();
        tmp += m * vertices[i].xyz;
        M += m;
    }

    tmp /= M;
    cm = tmp;
}

void Cell::setVisc(double mu)
{
    params.vertexVisc = mu / numberV;

    for (int i = 0; i < numberV; i++)
    {
        vertices[i].setVisc(params.vertexVisc);
    }

    params.totalVisc = getVisc();
}

double Cell::getVisc()
{
    double v = 0;

    for (int i = 0; i < numberV; i++)
    {
        v += vertices[i].getVisc();
    }

    return v;
}

double Cell::getMass()
{
    double totalMass = 0.0;

    for (int i = 0; i < numberV; i++)
    {
        totalMass += vertices[i].getMass();
    }

    return totalMass;
}

void Cell::setMass(double totm)
{
    params.vertexMass = totm / numberV;

    for (int i = 0; i < numberV; i++)
    {
        vertices[i].setMass(params.vertexMass);
    }

    params.totalMass = getMass();
}

void Cell::addVelocity(const Vector3D& nv)
{
    for (int i = 0; i < numberV; i++)
    {
        vertices[i].velocity += nv;
    }
}

void Cell::addXYZ(const Vector3D& nxyz)
{
    for (int i = 0; i < numberV; i++)
    {
        vertices[i].xyz += nxyz;
    }
}

//void Cell::moveToXYZ(const Vector3D& nxyz)
//{
//    for (int i = 0; i < numberV; i++)
//    {
//        vertices[i].xyz = nxyz;
//    }
//}

int Cell::numberOfTris()
{
    return numberT;
}

int Cell::numberOfVerts()
{
    return numberV;
}

void Cell::setVertexR(double rv)
{
    params.r_vertex = rv;
}

void Cell::setEcc(double a)
{
    params.ecc = a;
}

void Cell::setDp(double dP)
{
    params.dp = dP;
}

void Cell::setGamma(double g)
{
    params.gamma = g;
}

void Cell::setCellId(int ix)
{
    cellId = ix;
}

void Cell::setVerletR(double vr)
{
    params.verletR = vr;
}

void Cell::setInitR(double rinit)
{
    params.init_r = rinit;
}

void Cell::setVolumeC(double vc)
{
    params.vc = vc;
}

void Cell::setGrowthRate(double gr)
{
    params.growth_rate = gr;
}

void Cell::setBudDiameter(double bd)
{
    params.bud_d = bd;
}

void Cell::setDivisionRatio(double ds)
{
    params.div_ratio = ds;
}

void Cell::setNRT(double dp)
{
    nRT = dp * ( calcVolume() - OsmoticForce::getEpsilon() );
}

double Cell::getInitR()
{
    return params.init_r;
}

Vector3D Cell::getCm()
{
    calcCM();
    return cm;
}

double Cell::getVertexR()
{
    return params.r_vertex;
}

Vector3D& Cell::getVertexXYZ(int idx)
{
    return vertices[idx].xyz;
}

Vector3D& Cell::getVertexForce(int idx)
{
    return vertices[idx].force;
}

void Cell::getDistance(Vector3D& dkj, const Vector3D& vj, const Vector3D& vk, Box& box)
{
    dkj = vk - vj;

    if (box.pbc)
    {
        double x, y, z;
        double bsx = 2 * box.getX();
        double bsy = 2 * box.getY();
        double bsz = 2 * box.getZ();
        x = round(dkj.x / bsx) *  bsx;
        y = round(dkj.y / bsy) *  bsy;
        z = round(dkj.z / bsz) *  bsz;
        dkj.x -= x;
        dkj.y -= y;
        dkj.z -= z;
    }
}

void Cell::grow(double dt, double gr)
{  
    //std::cout << "numberV=" << numberV << std::endl;
    if (uniform() > 0.01)
        return;
        
    //double ptot = 0.0;
    //for (int i = 0; i < numberV; i++)
    //{
    //    ptot += (1.0 / vertices[i].numBonded);
    //}

    //double randn = uniform(0, 1.0);
    //double fracsum = 0.0;
    //int vertexId = -1;
    
    //for (int i = 0; i < numberV; i++)
    //{
    //    fracsum += (1. / vertices[i].numBonded) / (ptot);
    //    if (randn < fracsum)
    //    {
    //        vertexId = i;
    //        break;
    //    }
    //}
    
    int vertexId = getNextVertex();
    
    //vertexId = uniform(0, numberV);
    //vertexId = vcounter % numberV;
    //std::cout << "vertexId=" << vertexId << std::endl;
    //std::cout << "vertexId=" << vertexId << std::endl;
    //std::cout << " vertices[vertexId].numTris="<<vertices[vertexId].numTris << std::endl; 
    
    int triangle_num = uniform(0, vertices[vertexId].numTris);
    //std::cout << " triangle_num="<<triangle_num << std::endl; 
    
    
    int triangleId = vertices[vertexId].bondedTris[triangle_num];
    //std::cout << " triangle_id="<< triangleId << std::endl; 
    // check if the size of the triangle is large enough
    
    int vertPos = -1;
    int vert1 = -1;
    int vert2 = -1;
    if (triangles[triangleId].ia == vertexId)
    {
        vertPos = 0;
    }
    else if (triangles[triangleId].ib == vertexId)
    {
        vertPos = 1;
    }
    else
    {
         vertPos = 2;
    }
    
    if (vertPos == 0)
    {
        vert1 = triangles[triangleId].ib;
        vert2 = triangles[triangleId].ic;
    }
    else if(vertPos == 1)
    {
        vert1 = triangles[triangleId].ia;
        vert2 = triangles[triangleId].ic;
    }
    else
    {
        vert1 = triangles[triangleId].ia;
        vert2 = triangles[triangleId].ib;        
    }
    
    //std::cout << " vert1= "<< vert1 << " vert2=" << vert2 << std::endl;
    //std::cout << " vert1 tris#= "<< vertices[vert1].numTris << " vert2 tris3=" << vertices[vert2].numTris << std::endl;
    int secondTriangleId = -1;
    for (int i = 0; i < vertices[vert1].numTris; i++)
    {
        for (int j = 0; j < vertices[vert2].numTris; j++)
        {
            int t1 = vertices[vert1].bondedTris[i];
            int t2 = vertices[vert2].bondedTris[j];
            //std::cout <<  " t1="<<t1 <<" t2=" << t2 << std::endl;
            if (t1 == t2 && t1 != triangleId)
            {
                secondTriangleId = t1;
            }
        }
    }
    
    //std::cout << "secondTriangleId=" << secondTriangleId << std::endl;
    // ADD NEW VERTEX
    
    Vector3D newcoor = 0.5 * (vertices[vert1].xyz + vertices[vert2].xyz);
    vertices[numberV] = Vertex(newcoor.x, newcoor.y, newcoor.z);
    vertices[numberV].setId(numberV);
    vertices[numberV].setMass(vertices[vertexId].getMass());
    vertices[numberV].setVisc(vertices[vertexId].getVisc());
    int newid = numberV;
    numberV++;
    
    
    int vert3 = -1;
    if (triangles[secondTriangleId].ia == vert1)
    {
        if (triangles[secondTriangleId].ib == vert2)
        {
            vert3 = triangles[secondTriangleId].ic;
        } 
        else
        {
            vert3 = triangles[secondTriangleId].ib;
        }
    }

    if (triangles[secondTriangleId].ib == vert1)
    {
        if (triangles[secondTriangleId].ia == vert2)
        {
            vert3 = triangles[secondTriangleId].ic;
        } 
        else
        {
            vert3 = triangles[secondTriangleId].ia;
        }
    }    

    if (triangles[secondTriangleId].ic == vert1)
    {
        if (triangles[secondTriangleId].ia == vert2)
        {
            vert3 = triangles[secondTriangleId].ib;
        } 
        else
        {
            vert3 = triangles[secondTriangleId].ia;
        }
    }
    
    //triangles[triangleId].printVertexTriangle();
    //triangles[secondTriangleId].printVertexTriangle();
    
    //REMOVE OLD BONDS    
    vertices[vert1].removeNeighbor(vert2);
    vertices[vert2].removeNeighbor(vert1);
    
    // CREATE NEW BONDS
    
    Vector3D ab;// = vertices[aid].xyz - vertices[bid].xyz;
    ab = vertices[newid].xyz - vertices[vert1].xyz;
    //vertices[newid].addNeighbor(vert1, ab.length());
    //vertices[vert1].addNeighbor(newid, ab.length());
    
    ab = vertices[newid].xyz - vertices[vert2].xyz;
    //vertices[newid].addNeighbor(vert2, ab.length());
    //vertices[vert2].addNeighbor(newid, ab.length());
    
    ab = vertices[newid].xyz - vertices[vertexId].xyz;
    //vertices[newid].addNeighbor(vertexId, ab.length());
    //vertices[vertexId].addNeighbor(newid, ab.length());
    
    ab = vertices[newid].xyz - vertices[vert3].xyz;
    //vertices[newid].addNeighbor(vert3, ab.length() );
    //vertices[vert3].addNeighbor(newid, ab.length() );
    
    // =========
    ab = vertices[newid].xyz - vertices[vert1].xyz;
    vertices[newid].addNeighbor(vert1, ab.length() + 0.05 * r0av * uniform() );
    vertices[vert1].addNeighbor(newid, ab.length() + 0.05 * r0av * uniform() );
    
    ab = vertices[newid].xyz - vertices[vert2].xyz;
    vertices[newid].addNeighbor(vert2, ab.length() + 0.05*r0av * uniform() );
    vertices[vert2].addNeighbor(newid, ab.length() + 0.05*r0av * uniform() );
    
    ab = vertices[newid].xyz - vertices[vertexId].xyz;
    vertices[newid].addNeighbor(vertexId, ab.length() + 0.05*r0av * uniform() );
    vertices[vertexId].addNeighbor(newid, ab.length() + 0.05*r0av * uniform() );
    
    ab = vertices[newid].xyz - vertices[vert3].xyz;
    vertices[newid].addNeighbor(vert3, ab.length() + 0.05*r0av * uniform() );
    vertices[vert3].addNeighbor(newid, ab.length() + 0.05*r0av * uniform() );
    
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
    
    //std::cout << "trianlge 1 id=" << triangleId << " triangle 2 id=" <<secondTriangleId << std::endl;
   // std::cout << "vertexId=" << vertexId << " vert1=" << vert1 <<" vert2="<<vert2 << " vert3=" <<vert3 << " newid=" << newid<< std::endl;
    //std::cout << "r0av=" <<r0av << std::endl;
    
    // REMOVE TRIANGLES FROM THE VERTICES RECORD
    vertices[vert2].removeTriangle(triangleId);
    vertices[vert2].removeTriangle(secondTriangleId);
    
    triangles[triangleId].subsVertex(vert2, newid); 
    triangles[secondTriangleId].subsVertex(vert2, newid);
    vertices[newid].addTriangle(triangleId);
    vertices[newid].addTriangle(secondTriangleId);
    
    // NEW TRIANGLE 1
    VertexTriangle newTri1(vertexId, vert2, newid);
    triangles[numberT] = VertexTriangle(newTri1);
    triangles[numberT].setId(numberT);
    int tid1 = triangles[numberT].myindex;
    vertices[vertexId].addTriangle(tid1);
    vertices[vert2].addTriangle(tid1);
    vertices[newid].addTriangle(tid1);
    numberT++;
    
    
    // NEW TRIANGLE 2
    VertexTriangle newTri2(newid, vert2, vert3);
    triangles[numberT] = VertexTriangle(newTri2);
    triangles[numberT].setId(numberT);
    int tid2 = triangles[numberT].myindex;
    vertices[newid].addTriangle(tid2);
    vertices[vert2].addTriangle(tid2);
    vertices[vert3].addTriangle(tid2);
    numberT++;
    
    vcounter++;
    ///std::cout << "vcounter=" << vcounter << std::endl;


}

void Cell::growBud(double dt, double gr)
{
    
}

void Cell::divide()
{
    
}

void Cell::calcR0av()
{
    double r0_sum = 0.0;
    int N = 0;
    for (int i = 0; i < numberV; i++)
    {
        for (int j = 0; j < vertices[i].numBonded; j++)
        {
            r0_sum += vertices[i].r0[j];
            N++;
        }
    }
    r0av = r0_sum / N;
    for (int i = 0; i < numberV; i++)
    {
        //vertices[i].normalizedR0(r0av);
    }
    //std::cout << "bonds normalized" << std::endl;
}

void Cell::findBud()
{
    int round = 1;
    int triangles_added[MAX_T];
    int addedtris = 0;
    for (int i = 0; i < MAX_T; i++)
    {
        triangles_added[i] = -1;
    }
    // find and add a seed
    int budSeed = uniform(0, numberV);
    budIdx[budVno] = budSeed;
    budVno++;
    
    double bud_max_area = PI * params.bud_d * params.bud_d * 0.25;
    double bud_area = 0.0;
    
    bool uflag = false;
    
    while(true)
    {
       for (int i = 0; i < budVno; i++)
       {
           int vertid = budIdx[i];
           for (int j = 0; j < vertices[vertid].numTris; j++)
           {
               uflag = false;
               for (int k = 0; k < addedtris; k++)
               {
                   if (triangles_added[k] == vertices[vertid].bondedTris[j])
                   {
                       uflag = true;
                   }
               }
               if (!uflag)
               {
                   triangles_added[addedtris] = vertices[vertid].bondedTris[j];
                   addedtris++;
               }
           }
       }
       
       bud_area = 0.0;
       for (int i = 0; i < addedtris; i++)
       {
           int tri_ix = triangles_added[i];
           bud_area += triangles[tri_ix].area(vertices);
           std::cout <<" round #" << round << " tri idx=" << tri_ix << " tri area=" << triangles[tri_ix].area(vertices) << std::endl;
       }
       
       int ixa;
       int ixb;
       int ixc;
       
       bool flaga = false;
       bool flagb = false;
       bool flagc = false;
       
       if (bud_area <= bud_max_area)
       {
           // add unique verts
           for (int i = 0; i < addedtris; i++)
           {
               flaga = false;
               flagb = false;
               flagc = false;
               ixa = triangles[ triangles_added[i] ].ia;
               ixb = triangles[ triangles_added[i] ].ib;
               ixc = triangles[ triangles_added[i] ].ic;
               
               for (int j = 0; j < budVno; j++)
               {
                   if (ixa == budIdx[j])
                   {
                       flaga = true;
                   }
                   
                   if (ixb == budIdx[j])
                   {
                       flagb = true;
                   }
                   
                   if (ixc == budIdx[j])
                   {
                       flagc = true;
                   }
               }
               
               if (!flaga)
               {
                   budIdx[budVno] = ixa;
                   budVno++;
               }
               
               if (!flagb)
               {
                   budIdx[budVno] = ixb;
                   budVno++;
               }
               
               if (!flagc)
               {
                   budIdx[budVno] = ixc;
                   budVno++;
               }
           }
           round++;
       }
       else
       {
           break;
       }
    }
    
    std::cout << "bud_max_area=" << bud_max_area << " bud_area=" << bud_area << std::endl;
    for (int i = 0; i < budVno; i++)
    {
        std::cout << "i=" << i << " vno=" << budIdx[i] << std::endl;
    }
}

int Cell::getNextVertex()
{
    int vertexId = -1;
    // RANDOM
    return uniform(0, numberV);
    
    // 1 / #edges
    //double ptot = 0.0;
    //for (int i = 0; i < numberV; i++)
    //{
    //    ptot += (1.0 / vertices[i].numBonded);
    //}
    //double randn = uniform(0, 1.0);
    //double fracsum = 0.0;
    //
    //for (int i = 0; i < numberV; i++)
    //{
    //    fracsum += (1. / vertices[i].numBonded) / (ptot);
    //    if (randn < fracsum)
    //    {
    //        vertexId = i;
    //        break;
    //    }
    //}
    
    
    double spatialNb[numberV];
    for (int i = 0; i < numberV; i++)
    {
        spatialNb[i] = 0.0;
    }
    
    double ptot = 0.0;
    double cutoff = r0av;
    for (int i = 0; i< numberV; i++)
    {
        for (int j = 0; j < numberV; j++)
        {
            double r = (vertices[i].xyz - vertices[j].xyz).length() ;
            if (r <= cutoff && j!=i) 
            {
                spatialNb[i] += 1.0 / r;
                //spatialNb[j] += 1.0;// / r;
                //ptot += 1.0;// / r;
                //ptot += 1.0;// / r;
            }
        }
    }
    
    for (int i = 0; i < numberV; i++)
    {
        spatialNb[i] = std::max(spatialNb[i], 1.0);
    }
    
    for (int i = 0; i < numberV; i++)
    {
        ptot += 1.0 / spatialNb[i];
        std::cout << "spatialNb[" << i << "]=" << spatialNb[i] << std::endl;
    }
    std::cout << "ptot=" << ptot << std::endl;
    
    //for (int i = 0; i < numberV; i++)
    //{
    //    spatialNb[i] /= ptot;
    //}
    
    double sumcheck = 0.0;
    for (int i = 0; i < numberV; i++)
    {
        sumcheck += (1. / spatialNb[i]) / (ptot);
    }
    //std::cout << "sumcheck=" << sumcheck << std::endl;
    
    double randn = uniform(0, 1.0);
    double fracsum = 0.0;
   
    for (int i = 0; i < numberV; i++)
    {
        fracsum += (1. / spatialNb[i]) / (ptot);
        if (randn < fracsum)
        {
            vertexId = i;
            break;
        }
    }
    std::cout << "sumcheck=" << sumcheck <<  " fracsum="<<fracsum<< std::endl;
    std::cout << " vertexId=" << vertexId << std::endl;

    return vertexId;
}
