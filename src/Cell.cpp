#include "Cell.h"

utils::Logger Cell::cell_log("cell");

Cell::Cell(int depth) :  cellId(-1), my_phase(cell_phase_t::C_G1), numberV(0), numberT(0), nRT(0), r0av(0),
    vcounter(0), budVno(0)
{
    SimpleTriangulation sm(depth);
    std::list<Triangle> tris = sm.triangulate();
    Tinker::constructVertices(*this, tris);
    Tinker::constructVTriangles(*this, tris);
    Tinker::constructTopology(*this);
    calcAverageR0();
    randomRotate();
}

Cell::Cell(std::list<Triangle> tris) : cellId(-1), my_phase(cell_phase_t::C_G1), numberV(0), numberT(0), nRT(0),
    r0av(0), vcounter(0), budVno(0)
{
    Tinker::constructVertices(*this, tris);
    Tinker::constructVTriangles(*this, tris);
    Tinker::constructTopology(*this);
    calcAverageR0();
    randomRotate();
}

Cell::Cell(const Cell& orig) : cm_m(orig.cm_m), cm_b(orig.cm_b), vertices(orig.vertices), triangles(orig.triangles),
    cellId(orig.cellId), params(orig.params), my_phase(orig.my_phase), numberV(orig.numberV), numberT(orig.numberT), nRT(orig.nRT),
    r0av(orig.r0av), vcounter(orig.vcounter), budVno(orig.budVno)
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

double Cell::calcVolume(double eps)
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
        volume += tetr.volume(eps);
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

    if (Mb > 0.0)
    {
        tmp_b /= Mb;
        cm_b = tmp_b;
    }
}

void Cell::setVisc(double mu)
{
    params.vertexVisc = mu / numberV;
    //params.vertexVisc = mu;

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
    double correction = 4.0 * calcSurfaceArea() / sumL2();
    params.gamma = correction * g;
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
    setNRT(dp, 0.0);
    //nRT = dp * ( calcVolume() - OsmoticForce::getEpsilon() );
}

void Cell::setNRT(double dp, double ddp)
{
    double randu = uniform(0, ddp);
    nRT = (dp + randu) * ( calcVolume() - OsmoticForce::getEpsilon() );
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
    if (uniform() > params.growth_rate)
    {
        return;
    }

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

void Cell::calcAverageR0()
{
    double totSum = 0.0;

    for (int i = 0; i < numberV; i++)
    {
        for (int j = 0; j < vertices[i].numBonded; j++)
        {
            r0av += vertices[i].r0[j];
            totSum += 1.0;
        }
    }

    r0av /= totSum;
}

void Cell::setR0AvForAll()
{
    for (int i = 0; i < numberV; i++)
    {
        for (int j = 0; j < vertices[i].numBonded; j++)
        {
            vertices[i].r0[j] = r0av;
        }
    }
}

double Cell::sumL2()
{
    double sum_l2 = 0.0;

    for (int i = 0; i < numberV; i++)
    {
        for (int j = 0; j < vertices[i].numBonded; j++)
        {
            sum_l2 += vertices[i].r0[j] * vertices[i].r0[j];
        }
    }

    sum_l2 /= 2.0;
    return sum_l2;
}

double Cell::getPercLength(int i, int j)
{
    double r0 = vertices[i].r0[j];
    int k = vertices[i].bondedVerts[j];
    double r = (vertices[i].xyz - vertices[k].xyz).length();
    return 1.0 - r / r0;
}

double Cell::nbMagnitudeForce(std::vector<Cell> cells, Box& box, int vix)
{
    //std::cout << " vix=" << vix << std::endl;
    Vector3D force(0, 0, 0);

    for (int ci = 0; ci < cells.size(); ci++)
    {
        int ocellid = cells[ci].cellId;
        Vector3D dij;

        for (int i = 0; i < cells[ci].numberV; i++)
        {
            if (cellId != ocellid)
            {
                getDistance(dij, cells[ci].vertices[i].xyz, vertices[vix].xyz, box);
                //force += HertzianRepulsion::calcForce(dij, params.r_vertex, params.r_vertex, params.ecc);
            }
            else
            {
                if (i != vix && !vertices[i].isNeighbor(vix))
                {
                    getDistance(dij, cells[ci].vertices[i].xyz, vertices[vix].xyz, box);
                    //force += HertzianRepulsion::calcForce(dij, params.r_vertex, params.r_vertex, params.ecc);
                }
            }
        }
    }

    //Vector3D wallYZ, wallXZ, wallXY;
    //Vector3D dij;
    //double sgnx, sgny, sgnz;
    //double bsx = box.getX();
    //double bsy = box.getY();
    //double bsz = box.getZ();
    //double ecw = box.ecw;
    //for (int i = 0; i < numberV; i++)
    //{
    //sgnx = SIGN(vertices[vix].xyz.x);
    //wallYZ.x = sgnx * bsx;
    //wallYZ.y = vertices[vix].xyz.y;
    //wallYZ.z = vertices[vix].xyz.z;
    //dij = vertices[vix].xyz - wallYZ;
    //force += HertzianRepulsion::calcForce(dij, params.r_vertex, ecw);
    //sgny = SIGN(vertices[vix].xyz.y);
    //wallXZ.x = vertices[vix].xyz.x;
    //wallXZ.y = sgny * bsy;
    //wallXZ.z = vertices[vix].xyz.z;
    //dij = vertices[vix].xyz - wallXZ;
    //force += HertzianRepulsion::calcForce(dij, params.r_vertex, ecw);
    //sgnz = SIGN(vertices[vix].xyz.z);
    //wallXY.x = vertices[vix].xyz.x;
    //wallXY.y = vertices[vix].xyz.y;
    //wallXY.z = sgnz * bsz;
    //dij = vertices[vix].xyz - wallXY;
    //force += HertzianRepulsion::calcForce(dij, params.r_vertex, ecw);
    double vertexSurface = 0.0;

    for (int i = 0;  i < vertices[vix].numTris; i++)
    {
        vertexSurface += triangles[i].area(vertices);
        int j = vertices[vix].bondedTris[i];
        vertexSurface += triangles[j].area(vertices);
    }

    vertexSurface /= 3.0;
    //std::cout << "ar=" << 1.0 << std::endl;
    return force.length() / vertexSurface;
    //return force.length();
}

double Cell::nbMagnitudeForce(Cell ocell, Box& box)
{
    int ocellid = ocell.cellId;

    if (cellId == ocellid)
    {
        return 0.0;
    }

    Vector3D force(0, 0, 0);
    double totalContactSurface = 0.0;
    double total_force = 0.0;

    for (int vix = 0; vix < numberV; vix++)
    {
        force = Vector3D(0, 0, 0);
        Vector3D dij;

        for (int i = 0; i < ocell.numberV; i++)
        {
            getDistance(dij, ocell.vertices[i].xyz, vertices[vix].xyz, box);
            force += HertzianRepulsion::calcForce(dij, params.r_vertex, params.r_vertex, params.ecc);
        }

        total_force += force.length();

        if (force.length() > 0.0)
        {
            double vertexSurface = 0.0;

            for (int i = 0;  i < vertices[vix].numTris; i++)
            {
                vertexSurface += triangles[i].area(vertices);
                int j = vertices[vix].bondedTris[i];
                vertexSurface += triangles[j].area(vertices);
            }

            vertexSurface /= 3.0;
            totalContactSurface += vertexSurface;
        }
    }

    //return (total_force / totalContactSurface);
    return total_force;
}

void Cell::randomRotate()
{
    double u1 = uniform();
    double u2 = uniform();
    double u3 = uniform();

    double q0 = sqrt(1 - u1) * sin(2*M_PI*u2);
    double q1 = sqrt(1 - u1) * cos(2*M_PI*u2);
    double q2 = sqrt(u1) * sin(2*M_PI*u3);
    double q3 = sqrt(u1) * cos(2*M_PI*u3);

    double A[3][3];
    A[0][0] = q0*q0 + q1*q1 - q2*q2 - q3*q3;
    A[0][1] = 2*(q1*q2 + q0*q3);
    A[0][2] = 2*(q1*q3 - q0*q2);

    A[1][0] = 2*(q1*q2 - q0*q3);
    A[1][1] = q0*q0 - q1*q1 + q2*q2 - q3*q3;
    A[1][2] = 2*(q2*q3 + q0*q1);

    A[2][0] = 2*(q1*q3 + q0*q2);
    A[2][1] = 2*(q2*q3 - q0*q1);
    A[2][2] = q0*q0 - q1*q1 - q2*q2 + q3*q3;

    calcCM();

    double xnew = 0.0;
    double ynew = 0.0;
    double znew = 0.0;

    for (int i = 0; i  < numberV; i++)
    {
        double xi = vertices[i].xyz.x - cm_m.x;
        double yi = vertices[i].xyz.y - cm_m.y;
        double zi = vertices[i].xyz.z - cm_m.z;

        xnew  = A[0][0] * xi;
        xnew += A[0][1] * yi;
        xnew += A[0][2] * zi;

        ynew  = A[1][0] * xi;
        ynew += A[1][1] * yi;
        ynew += A[1][2] * zi;

        znew  = A[2][0] * xi;
        znew += A[2][1] * yi;
        znew += A[2][2] * zi;

        vertices[i].xyz.x = xnew + cm_m.x;
        vertices[i].xyz.y = ynew + cm_m.y;
        vertices[i].xyz.z = znew + cm_m.z;
    }
}

double Cell::contactForce(const Cell& other_cell, Box& box)
{
    int ocellid = other_cell.cellId;
    Vector3D dij;

    Vector3D force_collector(0, 0, 0);
    double contact_force = 0.0;

    for (int i = 0; i < numberV; i++)
    {
        for (int j = 0; j < other_cell.numberV; j++)
        {
            if (cellId != ocellid)
            {
                getDistance(dij, other_cell.vertices[j].xyz, vertices[i].xyz, box);
                force_collector += HertzianRepulsion::calcForce(dij, params.r_vertex, params.r_vertex, params.ecc);
            }
        }
        contact_force += force_collector.length();
        force_collector = Vector3D(0,0,0);
    }
    return contact_force;
}

//double Cell::contactForce(Box& box)
//{
//    Vector3D wallYZ, wallXZ, wallXY;
//    Vector3D dij;
//    double sgnx, sgny, sgnz;
//    double bsx = box.getX();
//    double bsy = box.getY();
//    double bsz = box.getZ();
//    double ecw = box.ecw;

//    Vector3D force_collector(0, 0, 0);
//    double contact_force = 0.0;

//    for (int i = 0; i < numberV; i++)
//    {
//        sgnx = SIGN(vertices[i].xyz.x);
//        wallYZ.x = sgnx * bsx;
//        wallYZ.y = vertices[i].xyz.y;
//        wallYZ.z = vertices[i].xyz.z;
//        dij = vertices[i].xyz - wallYZ;
//        force_collector += HertzianRepulsion::calcForce(dij, params.r_vertex, ecw);
//        sgny = SIGN(vertices[i].xyz.y);
//        wallXZ.x = vertices[i].xyz.x;
//        wallXZ.y = sgny * bsy;
//        wallXZ.z = vertices[i].xyz.z;
//        dij = vertices[i].xyz - wallXZ;
//        force_collector += HertzianRepulsion::calcForce(dij, params.r_vertex, ecw);
//        sgnz = SIGN(vertices[i].xyz.z);
//        wallXY.x = vertices[i].xyz.x;
//        wallXY.y = vertices[i].xyz.y;
//        wallXY.z = sgnz * bsz;
//        dij = vertices[i].xyz - wallXY;
//        force_collector += HertzianRepulsion::calcForce(dij, params.r_vertex, ecw);
//
//        contact_force += force_collector.length();
//        force_collector = Vector3D(0,0,0);
//    }
//    return force_collector;
//}

double Cell::contactArea(const Cell& ocell, Box& box)
{
    int ocellid = ocell.cellId;
    Vector3D dij;

    Vector3D force_collector1(0, 0, 0);
    Vector3D force_collector2(0, 0, 0);
    Vector3D force_collector3(0, 0, 0);

    int idx1, idx2, idx3;
    double fc1, fc2, fc3;

    double contact_area = 0.0;

    for (int i = 0; i < numberT; i++)
    {
        idx1 = triangles[i].ia;
        idx2 = triangles[i].ib;
        idx3 = triangles[i].ic;

        if (cellId != ocellid)
        {
            for (int j = 0; j < ocell.numberV; j++)
            {
                getDistance(dij, ocell.vertices[j].xyz, vertices[idx1].xyz, box);
                force_collector1 += HertzianRepulsion::calcForce(dij, params.r_vertex, params.r_vertex, params.ecc);

                getDistance(dij, ocell.vertices[j].xyz, vertices[idx2].xyz, box);
                force_collector2 += HertzianRepulsion::calcForce(dij, params.r_vertex, params.r_vertex, params.ecc);

                getDistance(dij, ocell.vertices[j].xyz, vertices[idx3].xyz, box);
                force_collector3 += HertzianRepulsion::calcForce(dij, params.r_vertex, params.r_vertex, params.ecc);
            }
        }

        fc1 = force_collector1.length();
        fc2 = force_collector2.length();
        fc3 = force_collector3.length();

        if (fc1*fc2*fc3 > 0)
        {
            contact_area += triangles[i].area(vertices);
        }

        force_collector1 = Vector3D(0,0,0);
        force_collector2 = Vector3D(0,0,0);
        force_collector3 = Vector3D(0,0,0);
        fc1 = 0;
        fc2 = 0;
        fc3 = 0;

    }

    return contact_area;
}

double Cell::contactArea(Box& box)
{
    Vector3D wallYZ, wallXZ, wallXY;
    Vector3D dij;
    double sgnx, sgny, sgnz;
    double bsx = box.getX();
    double bsy = box.getY();
    double bsz = box.getZ();
    double ecw = box.ecw;

    Vector3D force_collector(0, 0, 0);
    double contact_area = 0.0;

    int idxset [3] = {0, 0, 0};
    double fc [3] = {0, 0, 0};
    int idx;

    for (int i = 0; i < numberT; i++)
    {
        idxset[0] = triangles[i].ia;
        idxset[1] = triangles[i].ib;
        idxset[2] = triangles[i].ic;

        for (int j = 0; j < 3; j++)
        {
            idx = idxset[j];

            sgnx = SIGN(vertices[idx].xyz.x);
            wallYZ.x = sgnx * bsx;
            wallYZ.y = vertices[idx].xyz.y;
            wallYZ.z = vertices[idx].xyz.z;
            dij = vertices[idx].xyz - wallYZ;
            force_collector += HertzianRepulsion::calcForce(dij, params.r_vertex, ecw);

            sgny = SIGN(vertices[idx].xyz.y);
            wallXZ.x = vertices[idx].xyz.x;
            wallXZ.y = sgny * bsy;
            wallXZ.z = vertices[idx].xyz.z;
            dij = vertices[idx].xyz - wallXZ;
            force_collector += HertzianRepulsion::calcForce(dij, params.r_vertex, ecw);

            sgnz = SIGN(vertices[idx].xyz.z);
            wallXY.x = vertices[idx].xyz.x;
            wallXY.y = vertices[idx].xyz.y;
            wallXY.z = sgnz * bsz;
            dij = vertices[idx].xyz - wallXY;
            force_collector += HertzianRepulsion::calcForce(dij, params.r_vertex, ecw);

            fc[j] = force_collector.length();
            force_collector = Vector3D(0,0,0);
        }

        if (fc[0] * fc[1] * fc[2] > 0)
        {
            contact_area += triangles[i].area(vertices);
        }
        force_collector = Vector3D(0,0,0);
        fc[0] = 0;
        fc[1] = 0;
        fc[2] = 0;

    }
    return contact_area;
}

double Cell::surfaceStrainEnergy()
{
    double deps = 0.0;
    double R0ij;
    int idxj;

    for (int i = 0; i < numberV; i++)
    {
        for (int j = 0; j < vertices[i].numBonded; j++)
        {
            R0ij = vertices[i].r0[j];
            idxj = vertices[i].bondedVerts[j];
            //vertices[i].force += HookeanForce::calcForce(vertices[i].xyz, vertices[idxj].xyz, R0ij, params.gamma);
            Vector3D dR = vertices[idxj].xyz - vertices[i].xyz;
            double R = dR.length();
            deps += 0.5 * params.gamma * (R0ij - R) * (R0ij - R);
            //return f;
        }
    }
}