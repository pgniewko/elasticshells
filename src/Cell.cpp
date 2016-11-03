#include "Cell.h"

utils::Logger Cell::cell_log("cell");

bool Cell::no_bending = false;

Cell::Cell()
{

}

Cell::Cell(int depth)
{
    SimpleTriangulation sm(depth);
    std::list<Triangle> tris = sm.triangulate();
    Tinker::constructVertices(*this, tris);
    Tinker::constructVTriangles(*this, tris);
    Tinker::constructTopology(*this);
    Tinker::constructBSprings(*this);

}

Cell::Cell(std::list<Triangle> tris) : cell_id(-1), number_v(0), number_t(0), number_s(0), nRT(0),
    V0(0), fem_flag(false)//, bending_flag(true)
{
    Tinker::constructVertices(*this, tris);
    std::cout << "constructVertices DONE" << std::endl;
    Tinker::constructVTriangles(*this, tris);
    std::cout << "constructVTriangles DONE" << std::endl;
    Tinker::constructTopology(*this);
    std::cout << "constructTopology DONE" << std::endl;
    Tinker::constructBSprings(*this);
    std::cout << "constructBSprings DONE" << std::endl;
    randomRotate();

}

Cell::Cell(const Cell& orig) : cm_m(orig.cm_m), vertices(orig.vertices), triangles(orig.triangles), bhinges(orig.bhinges),
    cell_id(orig.cell_id), params(orig.params), number_v(orig.number_v), number_t(orig.number_t), number_s(orig.number_s),
    nRT(orig.nRT), V0(orig.V0),  fem_flag(orig.fem_flag) {} //, bending_flag(orig.bending_flag) {}

Cell::~Cell() {}

void Cell::calcBondedForces()
{
    if (fem_flag)
    {
        calcFemForces();
    }
    else
    {
        calcHarmonicForces();
    }

    calcOsmoticForces();
}

void Cell::calcHarmonicForces()
{
    double R0ij;
    double k0ij;
    int idxj;

    for (int i = 0; i < number_v; i++)
    {
        for (int j = 0; j < vertices[i].numBonded; j++)
        {
            R0ij = vertices[i].r0[j];
            k0ij = vertices[i].k0[j];
            idxj = vertices[i].bondedVerts[j];
            vertices[i].f_c += HookeanForce::calcForce(vertices[i].r_c, vertices[idxj].r_c, R0ij, k0ij);
        }
    }
}

void Cell::calcFemForces()
{

    for (int i = 0; i < number_t; i++)
    {
        triangles[i].calcFemForces(vertices);
    }

    if (!no_bending)
    {
        for (int i = 0; i < number_s; i++)
        {
            bhinges[i].calcBendingForces(vertices);
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
    double r1 = params.vertex_r;
    double r2 = other_cell.params.vertex_r;
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
                Box::getDistance(dij, other_cell.vertices[j].r_c, vertices[i].r_c, box);
                vertices[i].f_c += HertzianRepulsion::calcForce(dij, r1, r2, e1, e2, nu1, nu2);
            }
            else
            {
                if (i != j && !vertices[i].isNeighbor(j))
                {
                    Box::getDistance(dij, other_cell.vertices[j].r_c, vertices[i].r_c, box);
                    vertices[i].f_c += HertzianRepulsion::calcForce(dij, r1, r1, e1, e1, nu1, nu1);
                }
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
    double r1 = params.vertex_r;
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

double Cell::calcSurfaceArea(double d_param) const
{
    double totalSurface = 0.0;

    for (int i = 0; i < number_t; i++)
    {
        totalSurface += triangles[i].area(vertices, cm_m, d_param);
    }

    return totalSurface;
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
        volume += Tetrahedron::volume(vertices[va].r_c, vertices[vb].r_c, vertices[vc].r_c, cm_m, eps);
    }

    return volume;
}

void Cell::calcCM()
{
    Vector3D tmp_m(0.0, 0.0, 0.0);
    double Mm = 0.0;

    for (int i = 0; i < number_v; i++)
    {
        tmp_m += vertices[i].r_c;
        Mm += 1.0;
    }

    if (Mm > 0.0)
    {
        tmp_m /= Mm;
        cm_m = tmp_m;
    }
    else
    {
        // REPORT PROBLEM
    }
}

void Cell::setBSprings(double E, double t, double nu_)
{
    for (int i = 0; i < number_s; i++)
    {
        bhinges[i].setD(E, t, nu_);
        bhinges[i].setThetaZero(vertices);
    }
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
    params.vertex_r = rv;
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

void Cell::setSpringConst(double E, double t, double nu_, std::string model_t)
{
    if ( model_t.compare("ms_kot") == 0 )
    {
        double g = E * t * ( 2.0 / (1.0 - nu_) ) * calcSurfaceArea() / sumL2();

        for (int i = 0; i < number_v; i++ )
        {
            for (int j = 0; j < vertices[i].numBonded; j++)
            {
                vertices[i].k0[j] = g;
            }
        }
    }
    else if ( model_t.compare("ms_avg") == 0 )
    {
        int me, him, third;
        double area;

        int t_me;
        int t_him;

        double g_me_him;

        double a, b, c;

        for (int i = 0; i < number_v; i++)
        {
            me = i;

            for (int j = 0; j < vertices[me].numBonded; j++)
            {
                g_me_him = 0.0;
                him = vertices[me].bondedVerts[j];

                c = (vertices[me].r_c - vertices[him].r_c).length();

                for (int k = 0; k < vertices[me].numTris; k++)
                {
                    t_me = vertices[me].bondedTris[k];

                    for (int l = 0; l < vertices[him].numTris; l++)
                    {
                        t_him = vertices[him].bondedTris[l];

                        if (t_me == t_him)
                        {
                            third = -1;
                            area = triangles[t_me].area(vertices);

                            if (triangles[t_me].ia == me && triangles[t_me].ib == him)
                            {
                                third = triangles[t_me].ic;
                            }

                            if (triangles[t_me].ib == me && triangles[t_me].ia == him)
                            {
                                third = triangles[t_me].ic;
                            }

                            if (triangles[t_me].ia == me && triangles[t_me].ic == him)
                            {
                                third = triangles[t_me].ib;
                            }

                            if (triangles[t_me].ic == me && triangles[t_me].ia == him)
                            {
                                third = triangles[t_me].ib;
                            }

                            if (triangles[t_me].ib == me && triangles[t_me].ic == him)
                            {
                                third = triangles[t_me].ia;
                            }

                            if (triangles[t_me].ic == me && triangles[t_me].ib == him)
                            {
                                third = triangles[t_me].ia;
                            }

                            if (third == -1)
                            {
                                cell_log <<  utils::LogLevel::SEVERE  << "PROBLEM IN setSpringConst(). \n Simulation is exiting.\n";
                                exit(EXIT_FAILURE);
                            }


                            a = (vertices[me].r_c - vertices[third].r_c).length();
                            b = (vertices[third].r_c - vertices[him].r_c).length();

                            g_me_him += E * t * area / (c * c * (1 + nu_));
                            g_me_him += E * t * nu_ * (a * a + b * b - c * c) / ((1 - nu_ * nu_) * 8.0 * area);

                        }

                    }
                }

                vertices[me].k0[j] = g_me_him;
            }
        }
    }
    else if ( model_t.compare("fem") == 0 )
    {
        fem_flag = true;

        for (int i = 0; i < number_t; i++)
        {
            triangles[i].setParams(vertices, E, nu_, t);
        }
    }
    else
    {
        // print error and terminate
    }

}

void Cell::setCellId(int ix)
{
    cell_id = ix;

    for (int i = 0; i < number_v; i++)
    {
        vertices[i].setCellId(cell_id);
    }
}

void Cell::setInitR(double rinit)
{
    params.init_r = rinit;
}

double Cell::getInitR() const
{
    return params.init_r;
}

Vector3D Cell::getCm() const
{
    return cm_m;
}

double Cell::getVertexR() const
{
    return params.vertex_r;
}

double Cell::getE() const
{
    return params.ecc;
}

double Cell::getNu() const
{
    return params.nu;
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

    double u1 = uniform();
    double u2 = uniform();
    double u3 = uniform();
    double q0 = fastmath::fast_sqrt(1 - u1) * fastmath::fast_sin(2 * M_PI * u2);
    double q1 = fastmath::fast_sqrt(1 - u1) * fastmath::fast_cos(2 * M_PI * u2);
    double q2 = fastmath::fast_sqrt(u1) * fastmath::fast_sin(2 * M_PI * u3);
    double q3 = fastmath::fast_sqrt(u1) * fastmath::fast_cos(2 * M_PI * u3);
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
    Vector3D nj(0, 0, 0);
    int tj;

    for (int j = 0; j < vertices[vidx].numTris; j++)
    {
        tj = vertices[vidx].getTriangleId(j);

        if ( isInContact(tj, other_cell, box) )
        {
            nj = triangles[tj].normal(vertices);
            nj_fi = nj.x * force_collector.x + nj.y * force_collector.y + nj.z * force_collector.z;
            Aj = triangles[tj].area(vertices, cm_m, params.vertex_r);

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
            Aj = triangles[tj].area(vertices, cm_m, params.vertex_r);

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
    double r1 = params.vertex_r;
    double r2 = other_cell.params.vertex_r;
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
                Box::getDistance(dij, other_cell.vertices[j].r_c, vertices[i].r_c, box);
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

    double r1 = params.vertex_r;
    double r2 = other_cell.params.vertex_r;
    double e1 = params.ecc;
    double e2 = other_cell.params.ecc;
    double nu1 = params.nu;
    double nu2 = other_cell.params.nu;

    if (cell_id != ocellid)
    {
        for (int j = 0; j < other_cell.number_v; j++)
        {
            Box::getDistance(dij, other_cell.vertices[j].r_c, vertices[idx1].r_c, box);
            force_collector1 += HertzianRepulsion::calcForce(dij, r1, r2, e1, e2, nu1, nu2);

            Box::getDistance(dij, other_cell.vertices[j].r_c, vertices[idx2].r_c, box);
            force_collector2 += HertzianRepulsion::calcForce(dij, r1, r2, e1, e2, nu1, nu2);

            Box::getDistance(dij, other_cell.vertices[j].r_c, vertices[idx3].r_c, box);
            force_collector3 += HertzianRepulsion::calcForce(dij, r1, r2, e1, e2, nu1, nu2);
        }
    }

    fc1 = force_collector1.length_sq();
    fc2 = force_collector2.length_sq();
    fc3 = force_collector3.length_sq();

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

    Vector3D fc1v = box_force(box, idx1);
    Vector3D fc2v = box_force(box, idx2);
    Vector3D fc3v = box_force(box, idx3);

    double fc1 = fc1v.length_sq();
    double fc2 = fc2v.length_sq();
    double fc3 = fc3v.length_sq();
    
    if ( fc1*fc2*fc3 > 0)
    {
        return true;
    }

    return false;
}

double Cell::activeArea(const Box& box, const std::vector<Cell>& cells, double d_param) const
{
    double total_surface = calcSurfaceArea(d_param);

    double total_cell_cell_area = 0.0;
    for (uint cid = 0; cid < cells.size(); cid++)
    {
        if (cell_id != cells[cid].cell_id)
        {
            total_cell_cell_area += contactArea(cells[cid], box, d_param);
        }
    }

    double total_cell_box_area = contactArea2(box, d_param);
    double total_contact_area = total_cell_cell_area + total_cell_box_area;
    
    return std::max(0.0, total_surface - total_contact_area);
}

//double Cell::activeAreaFraction(const Box& box, const std::vector<Cell>& cells, double& counter, bool flag) const
//{
//    double total_surface = 0.0;
//
//    for (int t_idx = 0; t_idx < number_t; t_idx++)
//    {
//        total_surface += triangles[t_idx].area(vertices, cm_m, params.vertex_r);
//    }
//
//    double active_area = activeArea(box, cells, counter, flag);
//
//    return std::min(1.0, active_area / total_surface);
//}

double Cell::contactArea(const Cell& other_cell, const Box& box, const double d_param) const
{
    double dist = ( cm_m - other_cell.getCm() ).length();
    
    if (dist > 4.0 * params.init_r)
    {
        return 0.0;
    }
    
    double contact_area = 0.0;
    
    for (int t_idx = 0; t_idx < number_t; t_idx++)
    {
        if ( isInContact(t_idx, other_cell, box) )
        {
            contact_area += triangles[t_idx].area(vertices, cm_m, d_param);
        }
    }

    return contact_area;
}

double Cell::contactArea(const Cell& other_cell, const Box& box) const
{
    return contactArea(other_cell, box, params.vertex_r);
}

double Cell::contactArea(const Box& box, double d_param) const
{
    double contact_area = 0.0;

    for (int t_id = 0; t_id < number_t; t_id++)
    {
        if ( isInContact(t_id, box) )
        {
            // Two classes are affected by this code: CellBoxStress and WallCoverageFraction.
            double eps = params.vertex_r - d_param;

            if (eps < 0)
            {
                cell_log << utils::LogLevel::WARNING << "In contactArea(const Box& box, double d_param): ";
                cell_log << utils::LogLevel::WARNING << "d_param(=" << d_param << ") larger than params.vertex_r" << "\n";
            }

            contact_area += triangles[t_id].area(vertices, cm_m, eps);
        }
    }

    return contact_area;
}

double Cell::contactArea2(const Box& box, double d_param) const
{
    double contact_area = 0.0;

    for (int t_id = 0; t_id < number_t; t_id++)
    {
        if ( isInContact(t_id, box) )
        {
            // Two classes are affected by this code: ActiveActiveArea and ActiveActiveFraction.
            contact_area += triangles[t_id].area(vertices, cm_m, d_param);
        }
    }

    return contact_area;
}

double Cell::strainEnergy(const Box& box) const
{
    double deps = 0.0;
    double r0, r, k0;
    int neigh_idx;

    for (int i = 0; i < number_v; i++)
    {
        for (int j = 0; j < vertices[i].numBonded; j++)
        {
            r0 = vertices[i].r0[j];
            k0 = vertices[i].k0[j];
            neigh_idx = vertices[i].bondedVerts[j];
            Vector3D dR = vertices[neigh_idx].r_c - vertices[i].r_c;

            r = dR.length();
            deps += 0.5 * k0 * (r0 - r) * (r0 - r);
        }
    }

    return deps;
}

double Cell::maxStrain() const
{
    double max_strain = -10.0;
    double eps = 0.0;
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
            eps = (r - r0) / r0;

            max_strain = std::max(max_strain, eps );

        }
    }

    return max_strain;
}

double Cell::minStrain() const
{
    double max_strain = 10.0;
    double eps = 0.0;
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
            eps = (r - r0) / r0;

            max_strain = std::min(max_strain, eps);

        }
    }

    return max_strain;
}

double Cell::getStrain(int i, int j) const
{
    double r0 = vertices[i].r0[j];
    int neigh_idx = vertices[i].bondedVerts[j];
    Vector3D dR = vertices[neigh_idx].r_c - vertices[i].r_c;

    double r = dR.length();
    double eps = (r - r0) / r0;
    return eps;
}

double Cell::nbIntra(const Box& box) const
{
    return 0.0;
//    Vector3D dij;
//    Vector3D fij;
//    double r1 = params.vertex_r;
//    double e1 = params.ecc;
//    double nu1 = params.nu;
//
//    double nb_energy = 0.0;
//
//    int vertid;
//
//    for (int i = 0; i < number_v; i++)
//    {
//        for (int j = 0; j < vertices[i].numNbNeighs; j++)
//        {
//            if (vertices[i].nbCellsIdx[j] == cell_id)
//            {
//                vertid = vertices[i].nbVerts[j];
//                getDistance(dij, vertices[vertid].r_c, vertices[i].r_c, box);
//                fij = HertzianRepulsion::calcForce(dij, r1, r1, e1, e1, nu1, nu1);
//
//                if ( fij.length() > 0.0 )
//                {
//                    dij.set_length( 2 * r1 - dij.length() );
//                    nb_energy += fabs( fij.x * dij.x + fij.y * dij.y + fij.z * dij.z );
//                }
//            }
//        }
//    }
//
//    return nb_energy;
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

void Cell::setConstantVolume(double scale)
{
    params.vol_c = V0 * (scale * scale * scale);

    if (calcVolume() != V0)
    {
        cell_log << utils::LogLevel::WARNING << "(calcVolume() != V0) @setConstantVolume\n";
    }
}

double Cell::checkVolumeCondition(double eps)
{
    double V = calcVolume();
    return (params.vol_c - V) / V;
}

void Cell::ajustTurgor(double step)
{
    params.dp = (1.0 + step) * params.dp;
}

const cell_params_t& Cell::get_params() const
{
    return params;
}

std::ostream& operator<< (std::ostream& out, const Cell& c)
{
    out << "CELL " << c.cell_id << ' ';
    out << c.number_v << ' ' << c.number_t << ' ' << c.number_s << ' ';
    out << c.params.vertex_r << ' ' << c.params.ecc << ' ' << c.params.nu << ' ';
    out << c.params.dp << ' ' << c.params.init_r << ' ' << c.params.vol_c << ' ';
    out << c.nRT << ' ' << c.V0 << "\n";

    for (int i = 0; i < c.number_v; i++)
    {
        out << "CELLVERTEX " <<  c.cell_id << ' ' << c.vertices[i] << '\n';
    }

    for (int i = 0; i < c.number_t; i++)
    {
        out << "CELLTRIANG " <<  c.cell_id << ' ' << c.triangles[i] << '\n';
    }

    for (int i = 0; i < c.number_s; i++)
    {
        out << "CELLHINGE " <<  c.cell_id << ' ' << c.bhinges[i].getId() << ' ' << c.bhinges[i] << '\n';
    }



    return out;
}
