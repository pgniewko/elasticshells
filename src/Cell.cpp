#include "Cell.h"

utils::Logger Cell::cell_log("cell");

Cell::Cell(int depth) : cell_id(-1), my_phase(cell_phase_t::C_G1), number_v(0), number_t(0), nRT(0), V0(0),
    vert_no_bud(0)
{
    SimpleTriangulation sm(depth);
    std::list<Triangle> tris = sm.triangulate();
    Tinker::constructVertices(*this, tris);
    Tinker::constructVTriangles(*this, tris);
    Tinker::constructTopology(*this);
    randomRotate();
}

Cell::Cell(std::list<Triangle> tris) : cell_id(-1), my_phase(cell_phase_t::C_G1), number_v(0), number_t(0), nRT(0),
    V0(0), vert_no_bud(0)
{
    Tinker::constructVertices(*this, tris);
    Tinker::constructVTriangles(*this, tris);
    Tinker::constructTopology(*this);
    randomRotate();
}

Cell::Cell(const Cell& orig) : cm_m(orig.cm_m), cm_b(orig.cm_b), vertices(orig.vertices), triangles(orig.triangles),
    cell_id(orig.cell_id), params(orig.params), my_phase(orig.my_phase), number_v(orig.number_v), number_t(orig.number_t), nRT(orig.nRT),
    V0(orig.V0), vert_no_bud(orig.vert_no_bud)
{}

Cell::~Cell()
{}

void Cell::voidVerletLsit()
{
    for (int i = 0; i < number_v; i++)
    {
        vertices[i].numNbNeighs = 0;
    }
}

void Cell::builtVerletList(const Cell& other_cell, const Box& box)
{
    Vector3D distance_jk;
    double r_cut = 2 * params.r_vertex;

    if (this->cell_id != other_cell.cell_id)
    {
        for (int j = 0; j < number_v; j++)
        {
            for (int k = 0; k < other_cell.number_v; k++ )
            {
                getDistance(distance_jk,  other_cell.vertices[k].xyz, vertices[j].xyz, box);

                if (distance_jk.length() <= r_cut * params.verletR)
                {
                    vertices[j].addNbNeighbor(k, other_cell.cell_id);
                }
            }
        }
    }
    else
    {
        for (int j = 0; j < number_v; j++)
        {
            for (int k = 0; k < number_v; k++)
            {
                getDistance(distance_jk, vertices[k].xyz, vertices[j].xyz, box);

                if (j != k && !vertices[j].isNeighbor(k) && distance_jk.length() <= r_cut * params.verletR)
                {
                    vertices[j].addNbNeighbor(k, other_cell.cell_id);
                }
            }
        }
    }
}

void Cell::builtNbList(std::vector<Cell>& cells, DomainList& domains, const Box& box)
{
    int domainIdx;
    int domainn;
    Vector3D distance_ik;
    int vertIdx, cellIdx;
    double r_cut = 2 * params.r_vertex + EPSILON;

    for (int i = 0; i < number_v; i++)
    {
        domainIdx = domains.getDomainIndex(vertices[i]);

        for (int j = 0; j < domains.getNumberOfNeigh(domainIdx); j++)
        {
            domainn = domains.getDomainNeighbor(domainIdx, j);

            for (int k = 0; k < domains.getNumOfParticles(domainn); k++)
            {
                vertIdx = domains.getVertexIdx(domainn, k);
                cellIdx = domains.getCellIdx(domainn, k);

                if (this->cell_id != cellIdx)
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

    for (int i = 0; i < number_v; i++)
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

    for (int i = 0; i < number_v; i++)
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

    int iva, ivb, ivc;
    double turgor = getTurgor();

    for (int i = 0; i < number_t; i++)
    {
        iva = triangles[i].ia;
        ivb = triangles[i].ib;
        ivc = triangles[i].ic;
        Vector3D fa = OsmoticForce::calcForce(vertices[iva].xyz, vertices[ivb].xyz, vertices[ivc].xyz, cm_m, turgor);
        Vector3D fb = OsmoticForce::calcForce(vertices[ivb].xyz, vertices[ivc].xyz, vertices[iva].xyz, cm_m, turgor);
        Vector3D fc = OsmoticForce::calcForce(vertices[ivc].xyz, vertices[iva].xyz, vertices[ivb].xyz, cm_m, turgor);
        vertices[iva].force += fa;
        vertices[ivb].force += fb;
        vertices[ivc].force += fc;
    }
}

void Cell::calcNbForcesON2(const Cell& other_cell, const Box& box)
{
    int ocellid = other_cell.cell_id;
    Vector3D dij;
    double r1 = params.r_vertex;
    double r2 = other_cell.params.r_vertex;
    double e1 = params.ecc;
    double e2 = other_cell.params.ecc;
    double nu1 = params.nu;
    double nu2 = other_cell.params.nu;

    for (int i = 0; i < number_v; i++)
    {
        for (int j = 0; j < other_cell.number_v; j++)
        {
            if (cell_id != ocellid)
            {
                getDistance(dij, other_cell.vertices[j].xyz, vertices[i].xyz, box);
                vertices[i].force += HertzianRepulsion::calcForce(dij, r1, r2, e1, e2, nu1, nu2);
            }
            else
            {
                if (i != j && !vertices[i].isNeighbor(j))
                {
                    getDistance(dij, other_cell.vertices[j].xyz, vertices[i].xyz, box);
                    vertices[i].force += HertzianRepulsion::calcForce(dij, r1, r1, e1, e1, nu1, nu1);
                }
            }
        }
    }
}

void Cell::calcNbForcesVL(const Cell& other_cell, const Box& box)
{
    int ocellid = other_cell.cell_id;
    int vertid;
    Vector3D dij;
    double r1 = params.r_vertex;
    double r2 = other_cell.params.r_vertex;
    double e1 = params.ecc;
    double e2 = other_cell.params.ecc;
    double nu1 = params.nu;
    double nu2 = other_cell.params.nu;

    for (int i = 0; i < number_v; i++)
    {
        for (int j = 0; j < vertices[i].numNbNeighs; j++)
        {
            if (vertices[i].nbCellsIdx[j] == ocellid)
            {
                vertid = vertices[i].nbVerts[j];
                getDistance(dij, other_cell.vertices[vertid].xyz, vertices[i].xyz, box);
                vertices[i].force += HertzianRepulsion::calcForce(dij, r1, r2, e1, e2, nu1, nu2);
            }
        }
    }
}

void Cell::calcBoxForces(const Box& box)
{
    Vector3D wallYZ, wallXZ, wallXY;
    Vector3D dij;
    double sgnx, sgny, sgnz;
    double bsx = box.getX();
    double bsy = box.getY();
    double bsz = box.getZ();
    double eb  = box.getE();
    double nub = box.getNu();
    double rb_ = 0.0;
    double r1 = params.r_vertex;
    double e1 = params.ecc;
    double nu1 = params.nu;

    for (int i = 0; i < number_v; i++)
    {
        sgnx = SIGN(vertices[i].xyz.x);
        wallYZ.x = sgnx * bsx;
        wallYZ.y = vertices[i].xyz.y;
        wallYZ.z = vertices[i].xyz.z;
        dij = vertices[i].xyz - wallYZ;
        vertices[i].force += HertzianRepulsion::calcForce(dij, r1, rb_, e1, eb, nu1, nub);
        sgny = SIGN(vertices[i].xyz.y);
        wallXZ.x = vertices[i].xyz.x;
        wallXZ.y = sgny * bsy;
        wallXZ.z = vertices[i].xyz.z;
        dij = vertices[i].xyz - wallXZ;
        vertices[i].force += HertzianRepulsion::calcForce(dij, r1, rb_, e1, eb, nu1, nub);
        sgnz = SIGN(vertices[i].xyz.z);
        wallXY.x = vertices[i].xyz.x;
        wallXY.y = vertices[i].xyz.y;
        wallXY.z = sgnz * bsz;
        dij = vertices[i].xyz - wallXY;
        vertices[i].force += HertzianRepulsion::calcForce(dij, r1, rb_, e1, eb, nu1, nub);
    }
}

void Cell::voidForces()
{
    for (int i = 0; i < number_v; i++)
    {
        vertices[i].voidForce();
    }
}

double Cell::calcSurfaceArea()
{
    double surface = 0.0;

    for (int i = 0; i < number_t; i++)
    {
        surface += triangles[i].area(vertices);
    }

    return surface;
}

double Cell::calcVolume()
{
    return calcVolume(0.0);
}

double Cell::calcVolume(double eps)
{
    double volume = 0.0;
    int va, vb, vc;
    calcCM();

    for (int i = 0; i < number_t; i++)
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

    for (int i = 0; i < number_v; i++)
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
    params.vertexVisc = mu / number_v;

    for (int i = 0; i < number_v; i++)
    {
        vertices[i].setVisc(params.vertexVisc);
    }

    params.totalVisc = getCellViscosity();
}

double Cell::getCellViscosity()
{
    double v = 0;

    for (int i = 0; i < number_v; i++)
    {
        v += vertices[i].getVisc();
    }

    return v;
}

double Cell::getMass()
{
    double totalMass = 0.0;

    for (int i = 0; i < number_v; i++)
    {
        totalMass += vertices[i].getMass();
    }

    return totalMass;
}

void Cell::setMass(double totm)
{
    params.vertexMass = totm / number_v;

    for (int i = 0; i < number_v; i++)
    {
        vertices[i].setMass(params.vertexMass);
    }

    params.totalMass = getMass();
}

void Cell::addVelocity(const Vector3D& nv)
{
    for (int i = 0; i < number_v; i++)
    {
        vertices[i].velocity += nv;
    }
}

void Cell::addXYZ(const Vector3D& nxyz)
{
    for (int i = 0; i < number_v; i++)
    {
        vertices[i].xyz += nxyz;
    }
}

int Cell::getNumberTriangles()
{
    return number_t;
}

int Cell::getNumberVertices()
{
    return number_v;
}

void Cell::setVertexR(double rv)
{
    params.r_vertex = rv;
}

void Cell::setEcc(double a)
{
    params.ecc = a;
}

void Cell::setNu(double nu)
{
    params.nu = nu;
}

void Cell::setDp(double dP)
{
    setDp(dP, 0.0);
}

void Cell::setDp(double dP, double ddp)
{
    double randu = uniform(-ddp, ddp);
    params.dp = dP + randu;
    V0 = calcVolume();
    nRT = params.dp * V0 * ( 1.0 - OsmoticForce::getEpsilon() );
}

void Cell::setSpringConst(double g)
{
    double correction = 4.0 * calcSurfaceArea() / sumL2();
    params.gamma = correction * g;
}

void Cell::setCellId(int ix)
{
    cell_id = ix;
}

void Cell::setVerletR(double vr)
{
    params.verletR = vr;
}

void Cell::setInitR(double rinit)
{
    params.init_r = rinit;
}

void Cell::setBuddingVolume(double vc)
{
    params.div_volume = vc;
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

double Cell::getE()
{
    return params.ecc;
}

double Cell::getNu()
{
    return params.nu;
}

Vector3D& Cell::getVertexXYZ(int idx)
{
    return vertices[idx].xyz;
}

Vector3D& Cell::getVertexForce(int idx)
{
    return vertices[idx].force;
}

void Cell::getDistance(Vector3D& dkj, const Vector3D& vj, const Vector3D& vk, const Box& box)
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

double Cell::sumL2()
{
    double sum_l2 = 0.0;

    for (int i = 0; i < number_v; i++)
    {
        for (int j = 0; j < vertices[i].numBonded; j++)
        {
            sum_l2 += vertices[i].r0[j] * vertices[i].r0[j];
        }
    }

    sum_l2 /= 2.0;
    return sum_l2;
}

void Cell::randomRotate()
{
    double u1 = uniform();
    double u2 = uniform();
    double u3 = uniform();
    double q0 = sqrt(1 - u1) * sin(2 * M_PI * u2);
    double q1 = sqrt(1 - u1) * cos(2 * M_PI * u2);
    double q2 = sqrt(u1) * sin(2 * M_PI * u3);
    double q3 = sqrt(u1) * cos(2 * M_PI * u3);
    double A[3][3];
    A[0][0] = q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3;
    A[0][1] = 2 * (q1 * q2 + q0 * q3);
    A[0][2] = 2 * (q1 * q3 - q0 * q2);
    A[1][0] = 2 * (q1 * q2 - q0 * q3);
    A[1][1] = q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3;
    A[1][2] = 2 * (q2 * q3 + q0 * q1);
    A[2][0] = 2 * (q1 * q3 + q0 * q2);
    A[2][1] = 2 * (q2 * q3 - q0 * q1);
    A[2][2] = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;
    calcCM();
    double xnew = 0.0;
    double ynew = 0.0;
    double znew = 0.0;

    for (int i = 0; i  < number_v; i++)
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

//double Cell::contactForce(const Cell& other_cell, Box& box)
//{
//    int ocellid = other_cell.cell_id;
//    Vector3D dij;
//    Vector3D force_collector(0, 0, 0);
//    double contact_force = 0.0;
//    double r1 = params.r_vertex;
//    double r2 = other_cell.params.r_vertex;
//    double e1 = params.ecc;
//    double e2 = other_cell.params.ecc;
//    double nu1 = params.nu;
//    double nu2 = other_cell.params.nu;
//
//    for (int i = 0; i < number_v; i++)
//    {
//        for (int j = 0; j < other_cell.number_v; j++)
//        {
//            if (cell_id != ocellid)
//            {
//                getDistance(dij, other_cell.vertices[j].xyz, vertices[i].xyz, box);
//                force_collector += HertzianRepulsion::calcForce(dij, r1, r2, e1, e2, nu1, nu2);
//            }
//        }
//
//        contact_force += force_collector.length();
//        force_collector = Vector3D(0, 0, 0);
//    }
//
//    return contact_force;
//}

//double Cell::contactForceNew(const Cell& other_cell, Box& box)

double Cell::project_force(const Cell& other_cell, const Box& box, const Vector3D& force_collector, const int vidx)
{
    double fi = 0.0;
    double totAi = 0.0;
    double nj_fi = 0.0;
    double Aj = 0.0;
    Vector3D nj(0,0,0);
    int tj;

    for (int j = 0; j < vertices[vidx].numTris; j++)
    {
        tj = vertices[vidx].getTriangleId(j);
        if ( isInContact(tj, other_cell, box) )
        {
            nj = triangles[tj].normal(vertices);
            nj_fi = nj.x * force_collector.x + nj.y * force_collector.y + nj.z * force_collector.z;
            Aj = triangles[tj].area(vertices, cm_m, params.r_vertex);

            totAi += Aj;

            fi += fabs( nj_fi * Aj );
        }
    }

    if (totAi > 0)
    {
        fi /= totAi;
    }
    else
    {
        fi = 0.0;
    }
    
    return fi;
}

double Cell::project_force(const Box& box, const Vector3D& force_collector, const int vidx)
{
    double fi = 0.0;
    double totAi = 0.0;
    double nj_fi = 0.0;
    double Aj = 0.0;
    Vector3D nj;
    int tj;

    for (int j = 0; j < vertices[vidx].numTris; j++)
    {
        tj = vertices[vidx].getTriangleId(j);
        if ( isInContact(tj, box) )
        {
            nj = triangles[tj].normal(vertices);
            nj_fi = nj.x * force_collector.x + nj.y * force_collector.y + nj.z * force_collector.z;
            Aj = triangles[tj].area(vertices, cm_m, params.r_vertex);

            totAi += Aj;

            fi += fabs( nj_fi * Aj );
        }
    }

    if (totAi > 0)
    {
        fi /= totAi;
    }
    else
    {
        fi = 0.0;
    }
    
    return fi;
}

double Cell::contactForce(const Cell& other_cell, const Box& box)
{
    calcCM();
    int ocellid = other_cell.cell_id;
    Vector3D dij;
    Vector3D force_collector(0, 0, 0);
    double contact_force = 0.0;
    double r1 = params.r_vertex;
    double r2 = other_cell.params.r_vertex;
    double e1 = params.ecc;
    double e2 = other_cell.params.ecc;
    double nu1 = params.nu;
    double nu2 = other_cell.params.nu;

//    int tj;
    double fi;

//    Vector3D nj;
//    double nj_fi;
//    double totAi;
//    double Aj;

    for (int i = 0; i < number_v; i++)
    {
        for (int j = 0; j < other_cell.number_v; j++)
        {
            if (cell_id != ocellid)
            {
                getDistance(dij, other_cell.vertices[j].xyz, vertices[i].xyz, box);
                force_collector += HertzianRepulsion::calcForce(dij, r1, r2, e1, e2, nu1, nu2);
            }
        }
        
        fi = project_force(other_cell, box, force_collector, i);
        contact_force += fi;
        force_collector = Vector3D(0, 0, 0);

//        fi = 0.0;
//        totAi = 0.0;
//        nj_fi = 0.0;
//        Aj = 0.0;
//
//        // NEW CODE GOES HERE
//        for (int j = 0; j < vertices[i].numTris; j++)
//        {
//            tj = vertices[i].getTriangleId(j);
//
//            if ( isInContact(tj, other_cell, box) )
//            {
//                nj = triangles[tj].normal(vertices);
//                nj_fi = nj.x * force_collector.x + nj.y * force_collector.y + nj.z * force_collector.z;
//                Aj = triangles[tj].area(vertices, cm_m, params.r_vertex);
//
//                totAi += Aj;
//
//                fi += fabs( nj_fi * Aj );
//
//            }
//        }
//
//        if (totAi > 0)
//        {
//            fi /= totAi;
//        }
//        else
//        {
//            fi = 0.0;
//        }
    }

    return contact_force;
}

Vector3D Cell::box_force(const Box& box, const int vix)
{
    Vector3D wallYZ(0, 0, 0);
    Vector3D wallXZ(0, 0, 0);
    Vector3D wallXY(0, 0, 0);
    Vector3D force_collector(0, 0, 0);
    Vector3D djk(0, 0, 0);
    
    double sgnx, sgny, sgnz;
    double bsx = box.getX();
    double bsy = box.getY();
    double bsz = box.getZ();
    double eb = box.getE();
    double rb_ = 0.0;
    double nub = box.getNu();
    double e1 = getE();
    double r1 = getVertexR();
    double nu1 = getNu();
    
    Vector3D vertXYZ = getVertexXYZ(vix);
    
    sgnx = SIGN(vertXYZ.x);
    wallYZ.x = sgnx * bsx;
    wallYZ.y = vertXYZ.y;
    wallYZ.z = vertXYZ.z;
    djk = vertXYZ - wallYZ;
    force_collector += HertzianRepulsion::calcForce(djk, r1, rb_, e1, eb, nu1, nub);

    sgny = SIGN(vertXYZ.y);
    wallXZ.x = vertXYZ.x;
    wallXZ.y = sgny * bsy;
    wallXZ.z = vertXYZ.z;
    djk = vertXYZ - wallXZ;
    force_collector += HertzianRepulsion::calcForce(djk, r1, rb_, e1, eb, nu1, nub);

    sgnz = SIGN(vertXYZ.z);
    wallXY.x = vertXYZ.x;
    wallXY.y = vertXYZ.y;
    wallXY.z = sgnz * bsz;
    djk = vertXYZ - wallXY;
    force_collector += HertzianRepulsion::calcForce(djk, r1, rb_, e1, eb, nu1, nub);
        
    return force_collector;
}

double Cell::contactForce(const Box& box)
{
    if (box.pbc)
    {
        return 0.0;
    }

//    Vector3D wallYZ, wallXZ, wallXY;
//    Vector3D vertXYZ;
//    double sgnx, sgny, sgnz;
//    double bsx = box.getX();
//    double bsy = box.getY();
//    double bsz = box.getZ();
    //double fx, fy, fz;
    //Vector3D forceX(0, 0, 0);
    //Vector3D forceY(0, 0, 0);
    //Vector3D forceZ(0, 0, 0);

    Vector3D force_collector(0, 0, 0);
    //int tj;
    //double fi;

//    Vector3D nj;
    //double nj_fi;
    //double totAi;
    //double Aj;
    
    
    double contact_force = 0.0;
//    Vector3D djk;
//    double eb = box.getE();
//    double rb_ = 0.0;
//    double nub = box.getNu();
//    double e1 = getE();
//    double r1 = getVertexR();
//    double nu1 = getNu();


    //e1 = getE();
    //r1 = getVertexR();
    //nu1 = getNu();

    for (int i = 0; i < number_v; i++)
    {
        force_collector = box_force(box, i);
        contact_force += project_force(box, force_collector, i);
        
        //vertXYZ = getVertexXYZ(j);
        //sgnx = SIGN(vertXYZ.x);
        //wallYZ.x = sgnx * bsx;
        //wallYZ.y = vertXYZ.y;
        //wallYZ.z = vertXYZ.z;
        //djk = vertXYZ - wallYZ;
        //force_collector += HertzianRepulsion::calcForce(djk, r1, rb_, e1, eb, nu1, nub);
        //fx = forceX.length();
//
        //sgny = SIGN(vertXYZ.y);
        //wallXZ.x = vertXYZ.x;
        //wallXZ.y = sgny * bsy;
        //wallXZ.z = vertXYZ.z;
        //djk = vertXYZ - wallXZ;
        //force_collector += HertzianRepulsion::calcForce(djk, r1, rb_, e1, eb, nu1, nub);
        //fy = forceY.length();

        //sgnz = SIGN(vertXYZ.z);
        //wallXY.x = vertXYZ.x;
        //wallXY.y = vertXYZ.y;
        //wallXY.z = sgnz * bsz;
        //djk = vertXYZ - wallXY;
        //force_collector += HertzianRepulsion::calcForce(djk, r1, rb_, e1, eb, nu1, nub);
        //fz = forceZ.length();

        //totalForce +=  (fx + fy + fz);
        //contact_force += project_force(box, force_collector, j);
        
        //force_collector = Vector3D(0, 0, 0);
    }

    return contact_force;
}

double Cell::contactForceSF(const Box& box)
{
    if (box.pbc)
    {
        return 0.0;
    }

//    Vector3D wallYZ, wallXZ, wallXY;
//    Vector3D vertXYZ;
//    double sgnx, sgny, sgnz;
//    double bsx = box.getX();
//    double bsy = box.getY();
//    double bsz = box.getZ();
//    double fx, fy, fz;
//    Vector3D forceX(0, 0, 0);
//    Vector3D forceY(0, 0, 0);
//    Vector3D forceZ(0, 0, 0);
    Vector3D force_collector(0, 0, 0);
//    Vector3D djk;

    double contact_force = 0.0;
//    double eb = box.getE();
//    double rb_ = 0.0;
//    double nub = box.getNu();
//    double e1 = getE();
//    double r1 = getVertexR();
//    double nu1 = getNu();

    for (int i = 0; i < number_v; i++)
    {
        force_collector = box_force(box, i);
        contact_force += force_collector.length();
         
//        vertXYZ = getVertexXYZ(j);
//        sgnx = SIGN(vertXYZ.x);
//        wallYZ.x = sgnx * bsx;
//        wallYZ.y = vertXYZ.y;
//        wallYZ.z = vertXYZ.z;
//        djk = vertXYZ - wallYZ;
//        forceX = HertzianRepulsion::calcForce(djk, r1, rb_, e1, eb, nu1, nub);
//        fx = forceX.length();
//        sgny = SIGN(vertXYZ.y);
//        wallXZ.x = vertXYZ.x;
//        wallXZ.y = sgny * bsy;
//        wallXZ.z = vertXYZ.z;
//        djk = vertXYZ - wallXZ;
//        forceY = HertzianRepulsion::calcForce(djk, r1, rb_, e1, eb, nu1, nub);
//        fy = forceY.length();
//        sgnz = SIGN(vertXYZ.z);
//        wallXY.x = vertXYZ.x;
//        wallXY.y = vertXYZ.y;
//        wallXY.z = sgnz * bsz;
//        djk = vertXYZ - wallXY;
//        forceZ = HertzianRepulsion::calcForce(djk, r1, rb_, e1, eb, nu1, nub);
//        fz = forceZ.length();
//        contact_force += (fx + fy + fz);
    }

    return contact_force;
}

bool Cell::isInContact(int t_idx, const Cell& other_cell, const Box& box)
{
    int idx1, idx2, idx3;
    double fc1, fc2, fc3;

    int ocellid = other_cell.cell_id;

    Vector3D dij;
    Vector3D force_collector1(0, 0, 0);
    Vector3D force_collector2(0, 0, 0);
    Vector3D force_collector3(0, 0, 0);

    idx1 = triangles[t_idx].ia;
    idx2 = triangles[t_idx].ib;
    idx3 = triangles[t_idx].ic;

    double r1 = params.r_vertex;
    double r2 = other_cell.params.r_vertex;
    double e1 = params.ecc;
    double e2 = other_cell.params.ecc;
    double nu1 = params.nu;
    double nu2 = other_cell.params.nu;

    if (cell_id != ocellid)
    {
        for (int j = 0; j < other_cell.number_v; j++)
        {
            getDistance(dij, other_cell.vertices[j].xyz, vertices[idx1].xyz, box);
            force_collector1 += HertzianRepulsion::calcForce(dij, r1, r2, e1, e2, nu1, nu2);

            getDistance(dij, other_cell.vertices[j].xyz, vertices[idx2].xyz, box);
            force_collector2 += HertzianRepulsion::calcForce(dij, r1, r2, e1, e2, nu1, nu2);

            getDistance(dij, other_cell.vertices[j].xyz, vertices[idx3].xyz, box);
            force_collector3 += HertzianRepulsion::calcForce(dij, r1, r2, e1, e2, nu1, nu2);
        }
    }

    fc1 = force_collector1.length();
    fc2 = force_collector2.length();
    fc3 = force_collector3.length();

    if (fc1 * fc2 * fc3 > 0)
    {
        return true;
    }


    return false;
}

bool Cell::isInContact(int t_idx, const Box& box)
{
    if (box.pbc)
    {
        return false;
    }

    int idx1 = triangles[t_idx].ia;
    int idx2 = triangles[t_idx].ib;
    int idx3 = triangles[t_idx].ic;

    Vector3D fc1 = box_force(box, idx1);
    Vector3D fc2 = box_force(box, idx2);
    Vector3D fc3 = box_force(box, idx3);
    
    if ( (fc1.length() * fc2.length()* fc3.length()) > 0)
    {
        return true;
    }

    return false;
}

double Cell::contactArea(const Cell& other_cell, const Box& box)
{
    double contact_area = 0.0;
    calcCM();

    for (int i = 0; i < number_t; i++)
    {
        if ( isInContact(i, other_cell, box) )
        {
            contact_area += triangles[i].area(vertices, cm_m, params.r_vertex);
        }
    }

    return contact_area;
}

double Cell::contactArea(const Box& box, double d_param)
{
    calcCM();
    double contact_area = 0.0;
//    Vector3D wallYZ, wallXZ, wallXY;
//    Vector3D dij;
//    double sgnx, sgny, sgnz;
//    double bsx = box.getX();
//    double bsy = box.getY();
//    double bsz = box.getZ();
//    double eb  = box.getE();
//    double nub = box.getNu();
//    double rb_ = 0.0;
//    double r1 = params.r_vertex;
//    double e1 = params.ecc;
//    double nu1 = params.nu;
//    Vector3D force_collector(0, 0, 0);
//    double contact_area = 0.0;
//    int idxset [3] = {0, 0, 0};
//    double fc [3] = {0, 0, 0};
//    int idx;

    for (int t_idx = 0; t_idx < number_t; t_idx++)
    {
        if ( isInContact(t_idx, box) )
        {
            if (d_param > 0.0)
            {
                contact_area += triangles[t_idx].area(vertices);
            }
            else
            {
                contact_area += triangles[t_idx].area(vertices, cm_m, params.r_vertex);
            }            
        }
//        
//        idxset[0] = triangles[i].ia;
//        idxset[1] = triangles[i].ib;
//        idxset[2] = triangles[i].ic;
//
//        for (int j = 0; j < 3; j++)
//        {
//            idx = idxset[j];
//            sgnx = SIGN(vertices[idx].xyz.x);
//            wallYZ.x = sgnx * bsx;
//            wallYZ.y = vertices[idx].xyz.y;
//            wallYZ.z = vertices[idx].xyz.z;
//            dij = vertices[idx].xyz - wallYZ;
//            force_collector += HertzianRepulsion::calcForce(dij, r1, rb_, e1, eb, nu1, nub);
//            sgny = SIGN(vertices[idx].xyz.y);
//            wallXZ.x = vertices[idx].xyz.x;
//            wallXZ.y = sgny * bsy;
//            wallXZ.z = vertices[idx].xyz.z;
//            dij = vertices[idx].xyz - wallXZ;
//            force_collector += HertzianRepulsion::calcForce(dij, r1, rb_, e1, eb, nu1, nub);
//            sgnz = SIGN(vertices[idx].xyz.z);
//            wallXY.x = vertices[idx].xyz.x;
//            wallXY.y = vertices[idx].xyz.y;
//            wallXY.z = sgnz * bsz;
//            dij = vertices[idx].xyz - wallXY;
//            force_collector += HertzianRepulsion::calcForce(dij, r1, rb_, e1, eb, nu1, nub);
//            fc[j] = force_collector.length();
//            force_collector = Vector3D(0, 0, 0);
//        }
//
//        if (fc[0] * fc[1] * fc[2] > 0)
//        {
//            if (d_param > 0.0)
//            {
//                contact_area += triangles[i].area(vertices);
//            }
//            else
//            {
//                contact_area += triangles[i].area(vertices, cm_m, params.r_vertex);
//            }
//        }
//
//        force_collector = Vector3D(0, 0, 0);
//        fc[0] = 0;
//        fc[1] = 0;
//        fc[2] = 0;
    }

    return contact_area;
}

double Cell::surfaceStrainEnergy()
{
    double deps = 0.0;
    double R0ij, R;
    int idxj;

    for (int i = 0; i < number_v; i++)
    {
        for (int j = 0; j < vertices[i].numBonded; j++)
        {
            R0ij = vertices[i].r0[j];
            idxj = vertices[i].bondedVerts[j];
            Vector3D dR = vertices[idxj].xyz - vertices[i].xyz;
            R = dR.length();
            deps += 0.5 * params.gamma * (R0ij - R) * (R0ij - R);
        }
    }

    return deps;
}

double Cell::getTurgor()
{
    double turgor;

    if ( OsmoticForce::getFlag() )
    {
        double excluded_volume = V0 * OsmoticForce::getEpsilon();
        double cellVolume = calcVolume();

        if ( (cellVolume - excluded_volume) > 0 )
        {
            turgor = nRT / (cellVolume - excluded_volume);
        }
        else
        {
            cell_log << utils::LogLevel::SEVERE << " CELL VOLEUME SMALLER THAN OSMOTIC EXCLUDED VOLUME !\n";
            turgor = nRT / (0.01 * V0) ;
        }
    }
    else
    {
        turgor = params.dp;
    }

    return turgor;
}

// ********* CELL GROWTH
void Cell::cellCycle(double dt)
{
    switch (my_phase)
    {
        case cell_phase_t::C_G1:
            grow(dt);
            break;

        case cell_phase_t::C_SG2:
            bud(dt);
            break;

        case cell_phase_t::C_M:
            divide();
            break;

        default :
            break;

    }
}

void Cell::grow(double dt)
{
    if (uniform() > params.growth_rate * dt)
    {
        return;
    }

    Tinker::grow(*this);
}

void Cell::bud(double dt)
{
    findBud();

    if (uniform() > params.growth_rate * dt)
    {

        Tinker::bud(*this);
    }
}

void Cell::divide()
{
    Tinker::divide(*this);
}

void Cell::findBud()
{
    return;
}

// ********* END OF CELL GROWTH  --- DON'T ADD ANYTHIN BELOW