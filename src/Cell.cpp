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
                getDistance(distance_jk,  other_cell.vertices[k].r_c, vertices[j].r_c, box);

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
                getDistance(distance_jk, vertices[k].r_c, vertices[j].r_c, box);

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
                    getDistance(distance_ik, cells[cellIdx].vertices[vertIdx].r_c, vertices[i].r_c, box);

                    if (distance_ik.length() <= r_cut)
                    {
                        vertices[i].addNbNeighbor(vertIdx, cellIdx);
                    }
                }
                else
                {
                    getDistance(distance_ik, vertices[vertIdx].r_c, vertices[i].r_c, box);

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
            vertices[i].f_c += HookeanForce::calcForce(vertices[i].r_c, vertices[idxj].r_c, R0ij, params.gamma);
        }
    }
}

void Cell::calcOsmoticForces()
{
    int iva, ivb, ivc;
    double turgor = getTurgor();

    for (int i = 0; i < number_t; i++)
    {
        iva = triangles[i].ia;
        ivb = triangles[i].ib;
        ivc = triangles[i].ic;
        Vector3D fa = OsmoticForce::calcForce(vertices[iva].r_c, vertices[ivb].r_c, vertices[ivc].r_c, cm_m, turgor);
        Vector3D fb = OsmoticForce::calcForce(vertices[ivb].r_c, vertices[ivc].r_c, vertices[iva].r_c, cm_m, turgor);
        Vector3D fc = OsmoticForce::calcForce(vertices[ivc].r_c, vertices[iva].r_c, vertices[ivb].r_c, cm_m, turgor);
        vertices[iva].f_c += fa;
        vertices[ivb].f_c += fb;
        vertices[ivc].f_c += fc;
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
                getDistance(dij, other_cell.vertices[j].r_c, vertices[i].r_c, box);
                vertices[i].f_c += HertzianRepulsion::calcForce(dij, r1, r2, e1, e2, nu1, nu2);
            }
            else
            {
                if (i != j && !vertices[i].isNeighbor(j))
                {
                    getDistance(dij, other_cell.vertices[j].r_c, vertices[i].r_c, box);
                    vertices[i].f_c += HertzianRepulsion::calcForce(dij, r1, r1, e1, e1, nu1, nu1);
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
                getDistance(dij, other_cell.vertices[vertid].r_c, vertices[i].r_c, box);
                vertices[i].f_c += HertzianRepulsion::calcForce(dij, r1, r2, e1, e2, nu1, nu2);
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
        sgnx = SIGN(vertices[i].r_c.x);
        wallYZ.x = sgnx * bsx;
        wallYZ.y = vertices[i].r_c.y;
        wallYZ.z = vertices[i].r_c.z;
        dij = vertices[i].r_c - wallYZ;
        vertices[i].f_c += HertzianRepulsion::calcForce(dij, r1, rb_, e1, eb, nu1, nub);
        sgny = SIGN(vertices[i].r_c.y);
        wallXZ.x = vertices[i].r_c.x;
        wallXZ.y = sgny * bsy;
        wallXZ.z = vertices[i].r_c.z;
        dij = vertices[i].r_c - wallXZ;
        vertices[i].f_c += HertzianRepulsion::calcForce(dij, r1, rb_, e1, eb, nu1, nub);
        sgnz = SIGN(vertices[i].r_c.z);
        wallXY.x = vertices[i].r_c.x;
        wallXY.y = vertices[i].r_c.y;
        wallXY.z = sgnz * bsz;
        dij = vertices[i].r_c - wallXY;
        vertices[i].f_c += HertzianRepulsion::calcForce(dij, r1, rb_, e1, eb, nu1, nub);
    }
}

void Cell::voidForces()
{
    for (int i = 0; i < number_v; i++)
    {
        vertices[i].voidForce();
    }
}

double Cell::calcSurfaceArea() const
{
    double surface = 0.0;

    for (int i = 0; i < number_t; i++)
    {
        surface += triangles[i].area(vertices);
    }

    return surface;
}

double Cell::calcVolume(double eps) const
{
    double volume = 0.0;
    int va, vb, vc;

    for (int i = 0; i < number_t; i++)
    {
        va = triangles[i].ia;
        vb = triangles[i].ib;
        vc = triangles[i].ic;
        Tetrahedron tetr(vertices[va].r_c, vertices[vb].r_c, vertices[vc].r_c, cm_m);
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
            tmp_m += m * vertices[i].r_c;
            Mm += m;
        }

        if (vertices[i].getMyType() == vertex_t::BUD)
        {
            m = vertices[i].getMass();
            tmp_b += m * vertices[i].r_c;
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

double Cell::getCellViscosity() const
{
    double v = 0;

    for (int i = 0; i < number_v; i++)
    {
        v += vertices[i].getVisc();
    }

    return v;
}

void Cell::addXYZ(const Vector3D& nxyz)
{
    for (int i = 0; i < number_v; i++)
    {
        vertices[i].r_c += nxyz;
    }
}

int Cell::getNumberTriangles() const
{
    return number_t;
}

int Cell::getNumberVertices() const
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

double Cell::getInitR() const
{
    return params.init_r;
}

Vector3D Cell::getCm() const
{
//    calcCM();
    return cm_m;
}

double Cell::getVertexR() const
{
    return params.r_vertex;
}

double Cell::getE() const
{
    return params.ecc;
}

double Cell::getNu() const
{
    return params.nu;
}

void Cell::getDistance(Vector3D& dkj, const Vector3D& vj, const Vector3D& vk, const Box& box) const
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

double Cell::sumL2() const
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
    calcCM();
    //  
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
    double xnew = 0.0;
    double ynew = 0.0;
    double znew = 0.0;

    for (int i = 0; i  < number_v; i++)
    {
        double xi = vertices[i].r_c.x - cm_m.x;
        double yi = vertices[i].r_c.y - cm_m.y;
        double zi = vertices[i].r_c.z - cm_m.z;
        xnew  = A[0][0] * xi;
        xnew += A[0][1] * yi;
        xnew += A[0][2] * zi;
        ynew  = A[1][0] * xi;
        ynew += A[1][1] * yi;
        ynew += A[1][2] * zi;
        znew  = A[2][0] * xi;
        znew += A[2][1] * yi;
        znew += A[2][2] * zi;
        vertices[i].r_c.x = xnew + cm_m.x;
        vertices[i].r_c.y = ynew + cm_m.y;
        vertices[i].r_c.z = znew + cm_m.z;
    }
}

double Cell::project_force(const Cell& other_cell, const Box& box, const Vector3D& force_collector, const int vidx) const
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

double Cell::project_force(const Box& box, const Vector3D& force_collector, const int vidx) const
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

double Cell::contactForce(const Cell& other_cell, const Box& box) const
{
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
    
    double fi;

    for (int i = 0; i < number_v; i++)
    {
        for (int j = 0; j < other_cell.number_v; j++)
        {
            if (cell_id != ocellid)
            {
                getDistance(dij, other_cell.vertices[j].r_c, vertices[i].r_c, box);
                force_collector += HertzianRepulsion::calcForce(dij, r1, r2, e1, e2, nu1, nu2);
            }
        }
        
        fi = project_force(other_cell, box, force_collector, i);
        contact_force += fi;
        force_collector = Vector3D(0, 0, 0);
    }
    return contact_force;
}

Vector3D Cell::box_force(const Box& box, const int vix) const
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
    
    Vector3D vertXYZ = vertices[vix].r_c;
    
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

double Cell::contactForce(const Box& box) const
{
    if (box.pbc)
    {
        return 0.0;
    }

    Vector3D force_collector(0, 0, 0);
    double contact_force = 0.0;
    
    for (int i = 0; i < number_v; i++)
    {
        force_collector = box_force(box, i);
        contact_force += project_force(box, force_collector, i);
    }

    return contact_force;
}

double Cell::contactForceSF(const Box& box) const
{
    if (box.pbc)
    {
        return 0.0;
    }
    
    Vector3D force_collector(0, 0, 0);

    double contact_force = 0.0;

    for (int i = 0; i < number_v; i++)
    {
        force_collector = box_force(box, i);
        contact_force += force_collector.length();
    }

    return contact_force;
}

bool Cell::isInContact(int t_idx, const Cell& other_cell, const Box& box) const
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
            getDistance(dij, other_cell.vertices[j].r_c, vertices[idx1].r_c, box);
            force_collector1 += HertzianRepulsion::calcForce(dij, r1, r2, e1, e2, nu1, nu2);

            getDistance(dij, other_cell.vertices[j].r_c, vertices[idx2].r_c, box);
            force_collector2 += HertzianRepulsion::calcForce(dij, r1, r2, e1, e2, nu1, nu2);

            getDistance(dij, other_cell.vertices[j].r_c, vertices[idx3].r_c, box);
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

bool Cell::isInContact(int t_idx, const Box& box) const
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

double Cell::contactArea(const Cell& other_cell, const Box& box) const
{
    double contact_area = 0.0;

    for (int i = 0; i < number_t; i++)
    {
        if ( isInContact(i, other_cell, box) )
        {
            contact_area += triangles[i].area(vertices, cm_m, params.r_vertex);
        }
    }

    return contact_area;
}

double Cell::contactArea(const Box& box, double d_param) const
{
    double contact_area = 0.0;

    for (int t_idx = 0; t_idx < number_t; t_idx++)
    {
        if ( isInContact(t_idx, box) )
        {
            //TODO: this condition is confusing. 
            // It should be rather the other way around.
            // When changing be careful about changing observables.config parameters
            // Two classes are influenced by this code: CellBoxStress and WallCoverageFraction.
            if (d_param > 0.0)
            {
                contact_area += triangles[t_idx].area(vertices);
            }
            else
            {
                contact_area += triangles[t_idx].area(vertices, cm_m, params.r_vertex);
            }            
        }
    }

    return contact_area;
}

double Cell::strainEnergy(const Box& box) const
{
    double deps = 0.0;
    double r0, r;
    int neigh_idx;

    for (int i = 0; i < number_v; i++)
    {
        for (int j = 0; j < vertices[i].numBonded; j++)
        {
            r0 = vertices[i].r0[j];
            neigh_idx = vertices[i].bondedVerts[j];
            Vector3D dR = vertices[neigh_idx].r_c - vertices[i].r_c;
            
            r = dR.length();
            deps += 0.5 * params.gamma * (r0 - r) * (r0 - r);
        }
    }

    return deps;
}

double Cell::getTurgor() const
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

void Cell::update(double d)
{
    calcCM();
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