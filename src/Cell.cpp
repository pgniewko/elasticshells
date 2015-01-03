#include "Cell.h"

Cell::Cell(int depth) :  cellId(-1), my_phase(cell_phase_t::C_G0), numberV(0), numberT(0), nRT(0), r0av(0), 
        vcounter(0), budVno(0)
{
    SimpleTriangulation sm(depth);
    std::list<Triangle> tris = sm.triangulate();
    Tinker::constructVertices(*this, tris);
    Tinker::constructVTriangles(*this, tris);
    Tinker::constructTopology(*this);
}

Cell::Cell(std::list<Triangle> tris) : cellId(-1), my_phase(cell_phase_t::C_G0), numberV(0), numberT(0), nRT(0), 
        r0av(0), vcounter(0), budVno(0)
{
    Tinker::constructVertices(*this, tris);
    Tinker::constructVTriangles(*this, tris);
    Tinker::constructTopology(*this);
}

Cell::Cell(const Cell& orig) : cm_m(orig.cm_m), cm_b(orig.cm_b), vertices(orig.vertices), triangles(orig.triangles),
    cellId(orig.cellId), params(orig.params), my_phase(orig.my_phase), numberV(orig.numberV), numberT(orig.numberT), nRT(orig.nRT), r0av(orig.r0av),
    vcounter(orig.vcounter), budVno(orig.budVno)
{}

Cell::~Cell() {}

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
        Vector3D fa = OsmoticForce::calcForce(vertices[iva].xyz, vertices[ivb].xyz, vertices[ivc].xyz, cm_m, nRT, cellVolume, dp);
        Vector3D fb = OsmoticForce::calcForce(vertices[ivb].xyz, vertices[ivc].xyz, vertices[iva].xyz, cm_m, nRT, cellVolume, dp);
        Vector3D fc = OsmoticForce::calcForce(vertices[ivc].xyz, vertices[iva].xyz, vertices[ivb].xyz, cm_m, nRT, cellVolume, dp);
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
        Tetrahedron tetr(vertices[va].xyz, vertices[vb].xyz, vertices[vc].xyz, cm_m);
        volume += tetr.volume();
    }

    return volume;
}

void Cell::calcCM()
{
    Vector3D tmp_m(0.0, 0.0, 0.0);
    Vector3D tmp_b(0.0, 0.0, 0.0);
    double Mm = 0.0;
    double Mb = 0.0;
    double m;

    for (int i = 0; i < numberV; i++)
    {
        if (vertices[i].getMyType() == vertex_t::MOTHER)
        {
            m = vertices[i].getMass();
            tmp_m += m * vertices[i].xyz;
            Mm += m;
        }
        
        if (vertices[i].getMyType() == vertex_t::BUD)
        {
            m = vertices[i].getMass();
            tmp_b += m * vertices[i].xyz;
            Mb += m;
        }
    }

    if (Mm > 0.0)
    {
        tmp_m /= Mm;
        cm_m = tmp_m;
    }
    
    if (Mb > 0.0){
        tmp_b /= Mb;
        cm_b = tmp_b;
    }
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
    return cm_m;
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

void Cell::cellCycle()
{
    
    if (my_phase == cell_phase_t::C_G0)
    {
        // check if need to switch to C_G1
        // DO NOTHING
    }
    else if (my_phase == cell_phase_t::C_G1)
    {
        grow();
    }
    
    else if (my_phase == cell_phase_t::C_SG2)
    {
        bud();
    }
    else if (my_phase == cell_phase_t::C_M)
    {
        divide();
    }
}

void Cell::grow()
{
    if (uniform() > 0.00)
        return;
        
    Tinker::grow(*this);
}

void Cell::bud()
{
    findBud();
    Tinker::bud(*this);
}

void Cell::divide()
{
    Tinker::divide(*this);
}

void Cell::findBud()
{
    return;
}