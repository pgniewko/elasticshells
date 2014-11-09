#include "Cell.h"

Cell::Cell(int depth) :  cellId(-1), numberV(0), numberT(0), r_vertex(0),
    ecc(0), dp(0), gamma(0), verletR(0), initR(0), vertexVisc(0),
    vertexMass(0), totalVisc(0), totalMass(0), nRT(0)
{
    SimpleTriangulation sm(depth);
    std::list<Triangle> tris = sm.triangulate();
    constructVertices(tris);
    constructVTriangles(tris);
    constructTopology();
}

Cell::Cell(std::list<Triangle> tris) : cellId(-1), numberV(0), numberT(0), r_vertex(0),
    ecc(0), dp(0), gamma(0), verletR(0), initR(0), vertexVisc(0),
    vertexMass(0), totalVisc(0), totalMass(0),  nRT(0)
{
    constructVertices(tris);
    constructVTriangles(tris);
    constructTopology();
}

Cell::Cell(const Cell& orig) : cm(orig.cm), vertices(orig.vertices), triangles(orig.triangles),
    cellId(orig.cellId), numberV(orig.numberV), numberT(orig.numberT),
    r_vertex(orig.r_vertex), ecc(orig.ecc), dp(orig.dp), gamma(orig.gamma), verletR(orig.verletR),
    initR(orig.initR), vertexVisc(orig.vertexVisc), vertexMass(orig.vertexMass),
    totalVisc(orig.totalVisc), totalMass(orig.totalMass), nRT(orig.nRT)
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
    double r_cut = 2 * r_vertex;

    if (this->cellId != other_cell.cellId)
    {
        for (int j = 0; j < numberV; j++)
        {
            for (int k = 0; k < other_cell.numberV; k++ )
            {
                getDistance(distance_jk,  other_cell.vertices[k].xyz, vertices[j].xyz, box);

                if (distance_jk.length() <= r_cut * verletR)
                {
                    vertices[j].nbVerts[vertices[j].numNbNeighs] = k;
                    vertices[j].nbCellsIdx[vertices[j].numNbNeighs] = other_cell.cellId;
                    vertices[j].numNbNeighs++;
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

                if (j != k && !vertices[j].isNeighbor(k) && distance_jk.length() <= r_cut * verletR)
                {
                    vertices[j].nbVerts[vertices[j].numNbNeighs] = k;
                    vertices[j].nbCellsIdx[vertices[j].numNbNeighs] = other_cell.cellId;
                    vertices[j].numNbNeighs++;
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

    double r_cut = 2 * r_vertex;
    
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
                    { // TODO: REFACTOR adding nb vertex 
                        vertices[i].nbVerts[vertices[i].numNbNeighs] = vertIdx;
                        vertices[i].nbCellsIdx[vertices[i].numNbNeighs] = cellIdx;
                        vertices[i].numNbNeighs++;
                    }
                }
                else
                {
                    getDistance(distance_ik, vertices[vertIdx].xyz, vertices[i].xyz, box);

                    if (i != vertIdx && !vertices[i].isNeighbor(vertIdx) && distance_ik.length() <= r_cut)
                    {
                        vertices[i].nbVerts[vertices[i].numNbNeighs] = vertIdx;
                        vertices[i].nbCellsIdx[vertices[i].numNbNeighs] = cellIdx;
                        vertices[i].numNbNeighs++;
                    }
                }
            }
        }
    }
    
#ifdef TESTS
/* WRITE A CODE FOR SORTING NEIGHBORS*/
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
            vertices[i].force += HookeanForce::calcForce(vertices[i].xyz, vertices[idxj].xyz, R0ij, gamma);
        }
    }
}

void Cell::calcOsmoticForces()
{
    calcCM();
    double cellVolume = calcVolume();
    int iva, ivb, ivc;

    for (int i = 0; i < numberT; i++)
    {
        iva = triangles[i].ia;
        ivb = triangles[i].ib;
        ivc = triangles[i].ic;
        Vector3D fa = OsmoticForce::calcForce(vertices[iva].xyz, vertices[ivb].xyz, vertices[ivc].xyz, cm, nRT, cellVolume, dp);
        Vector3D fb = OsmoticForce::calcForce(vertices[ivb].xyz, vertices[ivc].xyz, vertices[iva].xyz, cm, nRT, cellVolume, dp);
        Vector3D fc = OsmoticForce::calcForce(vertices[ivc].xyz, vertices[iva].xyz, vertices[ivb].xyz, cm, nRT, cellVolume, dp);
        //Tetrahedron tetra(vertices[iva].xyz, vertices[ivb].xyz, vertices[ivc].xyz, cm);
        //Tetrahedron tetrb(vertices[ivb].xyz, vertices[ivc].xyz, vertices[iva].xyz, cm);
        //Tetrahedron tetrc(vertices[ivc].xyz, vertices[iva].xyz, vertices[ivb].xyz, cm);
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
                vertices[i].force += HertzianRepulsion::calcForce(dij, r_vertex, r_vertex, ecc);
            }
            else
            {
                if (i != j && !vertices[i].isNeighbor(j))
                {
                    getDistance(dij, other_cell.vertices[j].xyz, vertices[i].xyz, box);
                    vertices[i].force += HertzianRepulsion::calcForce(dij, r_vertex, r_vertex, ecc);
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
                vertices[i].force += HertzianRepulsion::calcForce(divix, r_vertex, r_vertex, ecc);
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
        vertices[i].force += HertzianRepulsion::calcForce(dij, r_vertex, ecw);
        
        sgny = SIGN(vertices[i].xyz.y);
        wallXZ.x = vertices[i].xyz.x;
        wallXZ.y = sgny * bsy;
        wallXZ.z = vertices[i].xyz.z;
        dij = vertices[i].xyz - wallXZ;
        vertices[i].force += HertzianRepulsion::calcForce(dij, r_vertex, ecw);
        
        sgnz = SIGN(vertices[i].xyz.z);
        wallXY.x = vertices[i].xyz.x;
        wallXY.y = vertices[i].xyz.y;
        wallXY.z = sgnz * bsz;
        dij = vertices[i].xyz - wallXY;
        vertices[i].force += HertzianRepulsion::calcForce(dij, r_vertex, ecw);
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
    vertexVisc = mu / numberV;

    for (int i = 0; i < numberV; i++)
    {
        vertices[i].setVisc(vertexVisc);
    }

    totalVisc = getVisc();
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
    vertexMass = totm / numberV;

    for (int i = 0; i < numberV; i++)
    {
        vertices[i].setMass(vertexMass);
    }

    totalMass = getMass();
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
    r_vertex = rv;
}

void Cell::setEcc(double a)
{
    ecc = a;
}

void Cell::setDp(double dP)
{
    dp = dP;
}

void Cell::setGamma(double g)
{
    gamma = g;
}

void Cell::setCellId(int ix)
{
    cellId = ix;
}

void Cell::setVerletR(double vr)
{
    verletR = vr;
}

void Cell::setInitR(double rinit)
{
    initR = rinit;
}

void Cell::setNRT(double dp)
{
    nRT = dp * ( calcVolume() - OsmoticForce::getEpsilon() );
}

double Cell::getInitR()
{
    return initR;
}

Vector3D Cell::getCm()
{
    calcCM();
    return cm;
}

double Cell::getVertexR()
{
    return r_vertex;
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