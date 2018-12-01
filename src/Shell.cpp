#include "Shell.h"

utils::Logger Shell::cell_log("shell");

bool Shell::no_bending = false;

double Shell::FORCE_FRAC(0.0);
double Shell::MIN_FORCE(0.0);
double Shell::MIN_FORCE_SQ(0.0);

Shell::Shell() {}

Shell::Shell(std::list<Triangle> tris) : shell_id(-1), number_v(0), number_t(0), number_s(0), nRT(0),
    V0(0)
{
    Tinker::constructVertices(*this, tris);
    Tinker::constructVTriangles(*this, tris);
    Tinker::constructTopology(*this);
    Tinker::constructBSprings(*this);
    randomRotate();
}

Shell::Shell(const Shell& orig) : center_of_mass(orig.center_of_mass), vertices(orig.vertices), triangles(orig.triangles), bhinges(orig.bhinges),
    shell_id(orig.shell_id), params(orig.params), number_v(orig.number_v), number_t(orig.number_t), number_s(orig.number_s),
    nRT(orig.nRT), V0(orig.V0) {}

Shell::~Shell() {}

//void Shell::calcBondedForces()
//{
//    if (number_v > 1 && number_t > 1)
//    {
//        calcFemForces();
//        calcOsmoticForces();
//    }
//}
//
//void Shell::calcFemForces()
//{
//
//    for (int i = 0; i < number_t; i++)
//    {
//        triangles[i].calcFemForces(vertices);
//    }
//
//    if (!no_bending)
//    {
//        for (int i = 0; i < number_s; i++)
//        {
//            bhinges[i].calcBendingForces(vertices);
//        }
//    }
//}
//
//void Shell::calcOsmoticForces()
//{
//    int iva, ivb, ivc;
//    double turgor = getTurgor();
//    
//    for (int i = 0; i < number_t; i++)
//    {
//        iva = triangles[i].ia;
//        ivb = triangles[i].ib;
//        ivc = triangles[i].ic;
//        Vector3D fa = OsmoticForce::calcForce(vertices[iva].r_c, vertices[ivb].r_c, vertices[ivc].r_c, center_of_mass, turgor);
//        Vector3D fb = OsmoticForce::calcForce(vertices[ivb].r_c, vertices[ivc].r_c, vertices[iva].r_c, center_of_mass, turgor);
//        Vector3D fc = OsmoticForce::calcForce(vertices[ivc].r_c, vertices[iva].r_c, vertices[ivb].r_c, center_of_mass, turgor);
//        vertices[iva].f_c += fa;
//        vertices[ivb].f_c += fb;
//        vertices[ivc].f_c += fc;
//    }
//}
//
//void Shell::calcNbForcesON2(const Shell& other_shell, const Box& box)
//{
//    int other_shell_id = other_shell.shell_id;
//    Vector3D dij;
//    double r1 = params.vertex_r;
//    double r2 = other_shell.params.vertex_r;
//    double e1 = params.ecc;
//    double e2 = other_shell.params.ecc;
//    double nu1 = params.nu;
//    double nu2 = other_shell.params.nu;
//
//    for (int i = 0; i < number_v; i++)
//    {
//        for (int j = 0; j < other_shell.number_v; j++)
//        {
//            if (shell_id != other_shell_id)
//            {
//                Box::getDistance(dij, other_shell.vertices[j].r_c, vertices[i].r_c, box);
//                vertices[i].f_c += HertzianRepulsion::calcForce(dij, r1, r2, e1, e2, nu1, nu2);
//            }
//            else
//            {
//                if (i != j && !vertices[i].isNeighbor(j))
//                {
//                    Box::getDistance(dij, other_shell.vertices[j].r_c, vertices[i].r_c, box);
//                    vertices[i].f_c += HertzianRepulsion::calcForce(dij, r1, r1, e1, e1, nu1, nu1);
//                }
//            }
//        }
//    }
//}
//
//void Shell::calcBoxForces(const Box& box)
//{
//    Vector3D wallYZ, wallXZ, wallXY;
//    Vector3D dij;
//    double sgnx, sgny, sgnz;
//    double bsx = box.getX();
//    double bsy = box.getY();
//    double bsz = box.getZ();
//    double eb  = box.getE();
//    double nub = box.getNu();
//    double rb_ = 0.0;
//    double r1 = params.vertex_r;
//    double e1 = params.ecc;
//    double nu1 = params.nu;
//
//    for (int i = 0; i < number_v; i++)
//    {
//        sgnx = SIGN(vertices[i].r_c.x);
//        wallYZ.x = sgnx * bsx;
//        wallYZ.y = vertices[i].r_c.y;
//        wallYZ.z = vertices[i].r_c.z;
//        dij = vertices[i].r_c - wallYZ;
//        vertices[i].f_c += HertzianRepulsion::calcForce(dij, r1, rb_, e1, eb, nu1, nub);
//        sgny = SIGN(vertices[i].r_c.y);
//        wallXZ.x = vertices[i].r_c.x;
//        wallXZ.y = sgny * bsy;
//        wallXZ.z = vertices[i].r_c.z;
//        dij = vertices[i].r_c - wallXZ;
//        vertices[i].f_c += HertzianRepulsion::calcForce(dij, r1, rb_, e1, eb, nu1, nub);
//        sgnz = SIGN(vertices[i].r_c.z);
//        wallXY.x = vertices[i].r_c.x;
//        wallXY.y = vertices[i].r_c.y;
//        wallXY.z = sgnz * bsz;
//        dij = vertices[i].r_c - wallXY;
//        vertices[i].f_c += HertzianRepulsion::calcForce(dij, r1, rb_, e1, eb, nu1, nub);
//    }
//}
//
//void Shell::voidForces()
//{
//    for (int i = 0; i < number_v; i++)
//    {
//        vertices[i].voidForce();
//    }
//}

double Shell::calcSurfaceArea(double d_param) const
{
    double totalSurface = 0.0;

    for (int i = 0; i < number_t; i++)
    {
        totalSurface += triangles[i].area(vertices, center_of_mass, d_param);
    }

    return totalSurface;
}

double Shell::calcSurfaceArea() const
{
    double surface = 0.0;

    for (int i = 0; i < number_t; i++)
    {
        surface += triangles[i].area(vertices);
    }

    return surface;
}

double Shell::calcVolume(double eps) const
{
    double volume = 0.0;

    if (number_v == 1)
    {
        volume = 4.0 / 3.0 * constants::pi * params.init_r * params.init_r * params.init_r;
    }
    else
    {
        int va, vb, vc;

        for (int i = 0; i < number_t; i++)
        {
            va = triangles[i].ia;
            vb = triangles[i].ib;
            vc = triangles[i].ic;
            volume += Tetrahedron::volume(vertices[va].r_c, vertices[vb].r_c, vertices[vc].r_c, center_of_mass, eps);
        }
    }

    return volume;
}

void Shell::calcCM()
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
        center_of_mass = tmp_m;
    }
    else
    {
        // REPORT PROBLEM
    }
}

void Shell::setBSprings(double E, double t, double nu_)
{
    if ( number_v == 1 || number_t == 1)
    {
        return;
    }

    for (int i = 0; i < number_s; i++)
    {
        bhinges[i].setD(E, t, nu_);
        bhinges[i].setThetaZero(vertices);
    }
}

void Shell::addXYZ(const Vector3D& nxyz)
{
    for (int i = 0; i < number_v; i++)
    {
        vertices[i].r_c += nxyz;
    }
}

int Shell::getNumberTriangles() const
{
    return number_t;
}

int Shell::getNumberVertices() const
{
    return number_v;
}

int Shell::getNumberHinges() const
{
    return number_s;
}

void Shell::setVertexR(double rv)
{
    params.vertex_r = rv;
}

void Shell::setEcc(double a)
{
    params.ecc = a;
}

void Shell::setNu(double nu)
{
    params.nu = nu;
}

void Shell::setDp(double dP)
{
    setDp(dP, 0.0);
}

void Shell::setDp(double dP, double ddp)
{
    double randu = uniform(-ddp, ddp);
    params.dp = dP + randu;
    V0 = calcVolume();
    nRT = params.dp * V0 * ( 1.0 - OsmoticForce::getEpsilon() );
}

void Shell::setSpringConst(double E, double t, double nu_, std::string model_t)
{
    for (int i = 0; i < number_t; i++)
    {
        triangles[i].setParams(vertices, E, nu_, t);
    }
}

void Shell::setShellId(int ix)
{
    shell_id = ix;

    for (int i = 0; i < number_v; i++)
    {
        vertices[i].set_shell_id(shell_id);
    }
}

void Shell::setInitR(double rinit)
{
    params.init_r = rinit;
}

double Shell::getInitR() const
{
    return params.init_r;
}

Vector3D Shell::getCm() const
{
    return center_of_mass;
}

double Shell::getVertexR() const
{
    return params.vertex_r;
}

double Shell::getE() const
{
    return params.ecc;
}

double Shell::getNu() const
{
    return params.nu;
}

void Shell::randomRotate()
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
        double xi = vertices[i].r_c.x - center_of_mass.x;
        double yi = vertices[i].r_c.y - center_of_mass.y;
        double zi = vertices[i].r_c.z - center_of_mass.z;
        xnew  = A[0][0] * xi;
        xnew += A[0][1] * yi;
        xnew += A[0][2] * zi;
        ynew  = A[1][0] * xi;
        ynew += A[1][1] * yi;
        ynew += A[1][2] * zi;
        znew  = A[2][0] * xi;
        znew += A[2][1] * yi;
        znew += A[2][2] * zi;
        vertices[i].r_c.x = xnew + center_of_mass.x;
        vertices[i].r_c.y = ynew + center_of_mass.y;
        vertices[i].r_c.z = znew + center_of_mass.z;
    }
}

//double Shell::project_force(const Shell& other_cell, const Box& box, const Vector3D& force_collector, const int vidx) const
//{
//    double fi = 0.0;
//    double totAi = 0.0;
//    double nj_fi = 0.0;
//    double Aj = 0.0;
//    Vector3D nj(0, 0, 0);
//    int tj;
//
//    for (int j = 0; j < vertices[vidx].facets_number; j++)
//    {
//        tj = vertices[vidx].getTriangleId(j);
//
//        if ( isInContact(tj, other_cell, box) )
//        {
//            nj = triangles[tj].normal(vertices);
//            nj_fi = nj.x * force_collector.x + nj.y * force_collector.y + nj.z * force_collector.z;
//            Aj = triangles[tj].area(vertices, center_of_mass, params.vertex_r);
//
//            totAi += Aj;
//
//            fi += fabs( nj_fi * Aj );
//        }
//    }
//
//    if (totAi > 0)
//    {
//        fi /= totAi;
//    }
//    else
//    {
//        fi = 0.0;
//    }
//
//    return fi;
//}
//
//double Shell::project_force(const Box& box, const Vector3D& force_collector, const int vidx) const
//{
//    double fi = 0.0;
//    double totAi = 0.0;
//    double nj_fi = 0.0;
//    double Aj = 0.0;
//    Vector3D nj;
//    int tj;
//
//    for (int j = 0; j < vertices[vidx].facets_number; j++)
//    {
//        tj = vertices[vidx].getTriangleId(j);
//
//        if ( isInContact(tj, box) )
//        {
//            nj = triangles[tj].normal(vertices);
//            nj_fi = nj.x * force_collector.x + nj.y * force_collector.y + nj.z * force_collector.z;
//            Aj = triangles[tj].area(vertices, center_of_mass, params.vertex_r);
//
//            totAi += Aj;
//
//            fi += fabs( nj_fi * Aj );
//        }
//    }
//
//    if (totAi > 0)
//    {
//        fi /= totAi;
//    }
//    else
//    {
//        fi = 0.0;
//    }
//
//    return fi;
//}
//
//Vector3D Shell::box_force(const Box& box, const int vix) const
//{
//    Vector3D wallYZ(0, 0, 0);
//    Vector3D wallXZ(0, 0, 0);
//    Vector3D wallXY(0, 0, 0);
//    Vector3D force_collector(0, 0, 0);
//    Vector3D djk(0, 0, 0);
//
//    double sgnx, sgny, sgnz;
//    double bsx = box.getX();
//    double bsy = box.getY();
//    double bsz = box.getZ();
//    double eb = box.getE();
//    double rb_ = 0.0;
//    double nub = box.getNu();
//    double e1 = getE();
//    double r1 = getVertexR();
//    double nu1 = getNu();
//
//    Vector3D vertXYZ = vertices[vix].r_c;
//
//    sgnx = SIGN(vertXYZ.x);
//    wallYZ.x = sgnx * bsx;
//    wallYZ.y = vertXYZ.y;
//    wallYZ.z = vertXYZ.z;
//    djk = vertXYZ - wallYZ;
//    force_collector += HertzianRepulsion::calcForce(djk, r1, rb_, e1, eb, nu1, nub);
//
//    sgny = SIGN(vertXYZ.y);
//    wallXZ.x = vertXYZ.x;
//    wallXZ.y = sgny * bsy;
//    wallXZ.z = vertXYZ.z;
//    djk = vertXYZ - wallXZ;
//    force_collector += HertzianRepulsion::calcForce(djk, r1, rb_, e1, eb, nu1, nub);
//
//    sgnz = SIGN(vertXYZ.z);
//    wallXY.x = vertXYZ.x;
//    wallXY.y = vertXYZ.y;
//    wallXY.z = sgnz * bsz;
//    djk = vertXYZ - wallXY;
//    force_collector += HertzianRepulsion::calcForce(djk, r1, rb_, e1, eb, nu1, nub);
//
//    return force_collector;
//}
//
//double Shell::contactForce(const Shell& other_cell, const Box& box, const bool flag) const
//{
//    int ocellid = other_cell.shell_id;
//    Vector3D dij;
//    Vector3D force_collector(0, 0, 0);
//    double contact_force = 0.0;
//    double r1 = params.vertex_r;
//    double r2 = other_cell.params.vertex_r;
//    double e1 = params.ecc;
//    double e2 = other_cell.params.ecc;
//    double nu1 = params.nu;
//    double nu2 = other_cell.params.nu;
//
//    double fi;
//
//    for (int i = 0; i < number_v; i++)
//    {
//        for (int j = 0; j < other_cell.number_v; j++)
//        {
//            if (shell_id != ocellid)
//            {
//                Box::getDistance(dij, other_cell.vertices[j].r_c, vertices[i].r_c, box);
//                force_collector += HertzianRepulsion::calcForce(dij, r1, r2, e1, e2, nu1, nu2);
//            }
//        }
//
//        if (flag)
//        {
//            contact_force += force_collector.length() ;
//        }
//        else
//        {
//            fi = project_force(other_cell, box, force_collector, i);
//            contact_force += fi;
//        }
//
//        force_collector = Vector3D(0, 0, 0);
//    }
//
//    return contact_force;
//}
//
//double Shell::contactForce(const Box& box) const
//{
//    if (box.pbc)
//    {
//        return 0.0;
//    }
//
//    Vector3D force_collector(0, 0, 0);
//    double contact_force = 0.0;
//
//    for (int i = 0; i < number_v; i++)
//    {
//        force_collector = box_force(box, i);
//        contact_force += project_force(box, force_collector, i);
//    }
//
//    return contact_force;
//}
//
//double Shell::contactForceSF(const Box& box) const // What class does it need ?
//{
//    if (box.pbc)
//    {
//        return 0.0;
//    }
//
//    Vector3D force_collector(0, 0, 0);
//
//    double contact_force = 0.0;
//
//    for (int i = 0; i < number_v; i++)
//    {
//        force_collector = box_force(box, i);
//        contact_force += force_collector.length();
//    }
//
//    return contact_force;
//}
//
//bool Shell::isInContact(int t_idx, const Shell& other_cell, const Box& box) const
//{
//    int idx1, idx2, idx3;
//    double fc1, fc2, fc3;
//
//    int ocellid = other_cell.shell_id;
//
//    Vector3D dij;
//    Vector3D force_collector1(0, 0, 0);
//    Vector3D force_collector2(0, 0, 0);
//    Vector3D force_collector3(0, 0, 0);
//
//    idx1 = triangles[t_idx].ia;
//    idx2 = triangles[t_idx].ib;
//    idx3 = triangles[t_idx].ic;
//
//    double r1 = params.vertex_r;
//    double r2 = other_cell.params.vertex_r;
//    double e1 = params.ecc;
//    double e2 = other_cell.params.ecc;
//    double nu1 = params.nu;
//    double nu2 = other_cell.params.nu;
//
//    if (shell_id != ocellid)
//    {
//        for (int j = 0; j < other_cell.number_v; j++)
//        {
//            Box::getDistance(dij, other_cell.vertices[j].r_c, vertices[idx1].r_c, box);
//            force_collector1 += HertzianRepulsion::calcForce(dij, r1, r2, e1, e2, nu1, nu2);
//
//            Box::getDistance(dij, other_cell.vertices[j].r_c, vertices[idx2].r_c, box);
//            force_collector2 += HertzianRepulsion::calcForce(dij, r1, r2, e1, e2, nu1, nu2);
//
//            Box::getDistance(dij, other_cell.vertices[j].r_c, vertices[idx3].r_c, box);
//            force_collector3 += HertzianRepulsion::calcForce(dij, r1, r2, e1, e2, nu1, nu2);
//        }
//    }
//
//    fc1 = force_collector1.length_sq();
//    fc2 = force_collector2.length_sq();
//    fc3 = force_collector3.length_sq();
//
//    if (fc1 > Shell::MIN_FORCE_SQ && fc2 > Shell::MIN_FORCE_SQ && fc3 > Shell::MIN_FORCE_SQ )
//    {
//        return true;
//    }
//
//    return false;
//}
//
//bool Shell::isInContact(int t_idx, const Box& box) const
//{
//    if (box.pbc)
//    {
//        return false;
//    }
//
//    int idx1 = triangles[t_idx].ia;
//    int idx2 = triangles[t_idx].ib;
//    int idx3 = triangles[t_idx].ic;
//
//    Vector3D fc1v = box_force(box, idx1);
//    Vector3D fc2v = box_force(box, idx2);
//    Vector3D fc3v = box_force(box, idx3);
//
//    double fc1 = fc1v.length_sq();
//    double fc2 = fc2v.length_sq();
//    double fc3 = fc3v.length_sq();
//
//    if (fc1 > Shell::MIN_FORCE_SQ && fc2 > Shell::MIN_FORCE_SQ && fc3 > Shell::MIN_FORCE_SQ )
//    {
//        return true;
//    }
//
//    return false;
//}
//
//bool Shell::isInContact(const Shell& other_cell, const Box& box) const
//{
//
//    if (shell_id == other_cell.shell_id)
//    {
//        return false;
//    }
//
//    Vector3D dij;
//    Vector3D force_ij(0, 0, 0);
//
//
//    double r1 = params.vertex_r;
//    double r2 = other_cell.params.vertex_r;
//    double e1 = params.ecc;
//    double e2 = other_cell.params.ecc;
//    double nu1 = params.nu;
//    double nu2 = other_cell.params.nu;
//
//
//    for (int i = 0; i < number_v; i++)
//    {
//        for (int j = 0; j < other_cell.number_v; j++)
//        {
//            Box::getDistance(dij, other_cell.vertices[j].r_c, vertices[i].r_c, box);
//            force_ij = HertzianRepulsion::calcForce(dij, r1, r2, e1, e2, nu1, nu2);
//
//            if (force_ij.length_sq() > Shell::MIN_FORCE_SQ)
//            {
//                return true;
//            }
//
//        }
//    }
//
//    return false;
//}
//
//double Shell::activeArea(const Box& box, const std::vector<Shell>& shells, double d_param) const
//{
//    double total_surface = calcSurfaceArea(d_param);
//
//    double total_cell_cell_area = 0.0;
//
//    for (uint cid = 0; cid < shells.size(); cid++)
//    {
//        if (shell_id != shells[cid].shell_id)
//        {
//            total_cell_cell_area += contactArea(shells[cid], box, d_param);
//        }
//    }
//
//    double total_cell_box_area = contactArea2(box, d_param);
//    double total_contact_area = total_cell_cell_area + total_cell_box_area;
//
//    return std::max(0.0, total_surface - total_contact_area);
//}
//
//double Shell::contactArea(const Shell& other_cell, const Box& box, const double d_param) const
//{
//    double dist = ( center_of_mass - other_cell.getCm() ).length();
//
//    if (dist > 4.0 * params.init_r)
//    {
//        return 0.0;
//    }
//
//    double contact_area = 0.0;
//
//    for (int t_idx = 0; t_idx < number_t; t_idx++)
//    {
//        if ( isInContact(t_idx, other_cell, box) )
//        {
//            contact_area += triangles[t_idx].area(vertices, center_of_mass, d_param);
//        }
//    }
//
//    return contact_area;
//}
//
//double Shell::contactArea(const Shell& other_cell, const Box& box) const
//{
//    return contactArea(other_cell, box, params.vertex_r);
//}
//
//double Shell::contactArea(const Box& box, double d_param) const
//{
//    double contact_area = 0.0;
//
//    for (int t_id = 0; t_id < number_t; t_id++)
//    {
//        if ( isInContact(t_id, box) )
//        {
//            // Two classes are affected by this code: CellBoxStress and WallCoverageFraction.
//            double eps = params.vertex_r - d_param;
//
//            if (eps < 0)
//            {
//                cell_log << utils::LogLevel::WARNING << "In contactArea(const Box& box, double d_param): ";
//                cell_log << utils::LogLevel::WARNING << "d_param(=" << d_param << ") larger than params.vertex_r" << "\n";
//            }
//
//            contact_area += triangles[t_id].area(vertices, center_of_mass, eps);
//        }
//    }
//
//    return contact_area;
//}
//
//double Shell::contactArea2(const Box& box, double d_param) const
//{
//    double contact_area = 0.0;
//
//    for (int t_id = 0; t_id < number_t; t_id++)
//    {
//        if ( isInContact(t_id, box) )
//        {
//            // Two classes are affected by this code: ActiveActiveArea and ActiveActiveFraction.
//            contact_area += triangles[t_id].area(vertices, center_of_mass, d_param);
//        }
//    }
//
//    return contact_area;
//}

double Shell::getTurgor() const
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

void Shell::update(double d)
{
    calcCM();
}

void Shell::setConstantVolume(double scale)
{
    params.vol_c = V0 * (scale * scale * scale);

    if (calcVolume() != V0)
    {
        cell_log << utils::LogLevel::WARNING << "(calcVolume() != V0) @setConstantVolume\n";
    }
}

double Shell::checkVolumeCondition()
{
    double V = calcVolume();
    return (params.vol_c - V) / V;
}

double Shell::ajust_turgor(double step)
{
    params.dp = (1.0 + step) * params.dp;
    return params.dp;
}

const shell_params_t& Shell::get_params() const
{
    return params;
}

std::ostream& operator<< (std::ostream& out, const Shell& c)
{
    out << "SHELL " << c.shell_id << ' ';
    out << c.number_v << ' ' << c.number_t << ' ' << c.number_s << ' ';
    out << c.params.vertex_r << ' ' << c.params.ecc << ' ' << c.params.nu << ' ';
    out << c.params.dp << ' ' << c.params.init_r << ' ' << c.params.vol_c << ' ';
    out << c.nRT << ' ' << c.V0 << "\n";

    for (int i = 0; i < c.number_v; i++)
    {
        out << "SHELLVERTEX " <<  c.shell_id << ' ' << c.vertices[i] << '\n';
    }

    for (int i = 0; i < c.number_t; i++)
    {
        out << "SHELLTRIANG " <<  c.shell_id << ' ' << c.triangles[i] << '\n';
    }

    for (int i = 0; i < c.number_s; i++)
    {
        out << "SHELLHINGE " <<  c.shell_id << ' ' << c.bhinges[i].getId() << ' ' << c.bhinges[i] << '\n';
    }

    return out;
}
