#include "ForcesCalculator.h"

ForcesCalculator::ForcesCalculator(int m_, bool pbc_, bool bend) : m(m_), pbc(pbc_), bending(bend), dl(m_, pbc_) 
{
    
}

ForcesCalculator::ForcesCalculator(const ForcesCalculator& orig) : m(orig.m), pbc(orig.pbc), bending(orig.bending), dl(orig.m, orig.pbc)
{
}

ForcesCalculator::~ForcesCalculator() {
}

void ForcesCalculator::calculate_forces(const std::vector<double>& xyz, 
                                        std::vector<double>& forces,
                                        const std::vector<element>& elements, 
                                        const std::vector<hinge>& hinges, 
                                        const std::vector<object_map> vs_map,
                                        const std::vector<double> turgors,
                                        const int num_shells) const
{
    
    
    // ITERATE OVER ELEMENTS
    
    // ITERATE OVER HINGES
    
    // CALCULATE MASS CENTERS
    // PRESSURE FORCES
    
    // END WITH NON-BONDED
    
    
}

void ForcesCalculator::evaluate_elements(const std::vector<double>& xyz,
                                         std::vector<double>& forces,
                                         const std::vector<element>& elements) const
{
    double x1, y1, z1;
    double x2, y2, z2;
    double x3, y3, z3;
    
    int vert_a, vert_b, vert_c;
    Vector3D va, vb, vc;
    
    double l0_sq, l1_sq, l2_sq;
    
    Vector3D T11, T12;
    Vector3D T21, T22;
    Vector3D T31, T32;
    
    element el;
    for (uint i = 0; i < elements.size(); i++)
    {
        el = elements[i];
        vert_a = el.ia;
        vert_b = el.ib;
        vert_c = el.ic;
        
        x1 = xyz[3 * vert_a + 0];
        y1 = xyz[3 * vert_a + 1];
        z1 = xyz[3 * vert_a + 2];
        
        x2 = xyz[3 * vert_b + 0];
        y2 = xyz[3 * vert_b + 1];
        z2 = xyz[3 * vert_b + 2];
        
        x3 = xyz[3 * vert_c + 0];
        y3 = xyz[3 * vert_c + 1];
        z3 = xyz[3 * vert_c + 2];
        
        va = Vector3D(x1, y1, z1);
        vb = Vector3D(x2, y2, z2);
        vc = Vector3D(x3, y3, z3);
        
        l0_sq = (vb - vc).length_sq() - el.L2[0];
        l1_sq = (va - vc).length_sq() - el.L2[1];
        l2_sq = (va - vb).length_sq() - el.L2[2];
        
        T11  =  el.ki[2] * l2_sq * (vb - va) + el.ki[1] * l1_sq * (vc - va);
        T12  = (el.ci[1] * l0_sq + el.ci[0] * l1_sq) * (vb - va);
        T12 += (el.ci[2] * l0_sq + el.ci[0] * l2_sq) * (vc - va);
        forces[3 * vert_a + 0] += (T11.x + T12.x);
        forces[3 * vert_a + 1] += (T11.y + T12.y);
        forces[3 * vert_a + 2] += (T11.z + T12.z);
        
        
        T21  =  el.ki[2] * l2_sq * (va - vb) + el.ki[0] * l0_sq * (vc - vb);
        T22  = (el.ci[0] * l1_sq + el.ci[1] * l0_sq) * (va - vb);
        T22 += (el.ci[2] * l1_sq + el.ci[1] * l2_sq) * (vc - vb);
        forces[3 * vert_b + 0] += (T21.x + T22.x);
        forces[3 * vert_b + 1] += (T21.y + T22.y);
        forces[3 * vert_b + 2] += (T21.z + T22.z);
 
        
        T31  = el.ki[1] * l1_sq * (va - vc) + el.ki[0] * l0_sq * (vb - vc);
        T32  = (el.ci[0] * l2_sq + el.ci[2] * l0_sq) * (va - vc);
        T32 += (el.ci[1] * l2_sq + el.ci[2] * l1_sq) * (vb - vc);        
        forces[3 * vert_c + 0] += (T31.x + T32.x);
        forces[3 * vert_c + 1] += (T31.y + T32.y);
        forces[3 * vert_c + 2] += (T31.z + T32.z);
    }
}


void ForcesCalculator::evaluate_hinges(const std::vector<double>& xyz,
                                       std::vector<double>& forces,
                                       const std::vector<hinge>& hinges) const
{
    double x1, y1, z1;
    double x2, y2, z2;
    double x3, y3, z3;
    double x4, y4, z4;
    
    int vert_1, vert_2, vert_3, vert_4;
    
    Vector3D ve1, ve2, ve3, ve4;
    Vector3D edge;
    Vector3D A1, A2;
    Vector3D n1, n2;
    
    Vector3D u1, u2, u3, u4;
    
    double edge_norm, edge_norm2;
    double area1, area2;
    double D, C, theta, theta0;
    hinge hi;
    for (uint i = 0; i < hinges.size(); i++)
    {
        hi = hinges[i];
        vert_1 = hi.v1;
        vert_2 = hi.v2;
        vert_3 = hi.v3;
        vert_4 = hi.v4;
        D = hi.D;
        theta0 = hi.theta0;
        
        x1 = xyz[3 * vert_1 + 0];
        y1 = xyz[3 * vert_1 + 1];
        z1 = xyz[3 * vert_1 + 2];
        
        x2 = xyz[3 * vert_2 + 0];
        y2 = xyz[3 * vert_2 + 1];
        z2 = xyz[3 * vert_2 + 2];
        
        x3 = xyz[3 * vert_3 + 0];
        y3 = xyz[3 * vert_3 + 1];
        z3 = xyz[3 * vert_3 + 2];
        
        x4 = xyz[3 * vert_4 + 0];
        y4 = xyz[3 * vert_4 + 1];
        z4 = xyz[3 * vert_4 + 2];
        
        ve1 = Vector3D(x1, y1, z1);
        ve2 = Vector3D(x2, y2, z2);
        ve3 = Vector3D(x3, y3, z3);
        ve4 = Vector3D(x4, y4, z4);
        
        edge = ve4 - ve3;
        edge_norm = edge.length();
        edge_norm2 = edge_norm * edge_norm;
        
        A1 = 0.5 * cross(ve1 - ve3, ve1 - ve4);
        area1 = A1.length();
        n1 = A1 / area1;

        A2 = 0.5 * cross(ve2 - ve4, ve2 - ve3);
        area2 = A2.length();
        n2 = A2 / area2;
        
        theta = calculate_theta(ve1, ve2, ve3, ve4);
        C = D * edge_norm2 / (area1 + area2) * fastmath::fast_sin(theta - theta0); // Grinspun eq.3, Cubic Shells, 2007
        
        u1 = edge_norm / (2.0 * area1) * n1;
        u2 = edge_norm / (2.0 * area2) * n2;

        u3 =  dot(ve1 - ve4, edge) * n1 / (2.0 * area1 * edge_norm) + dot(ve2 - ve4, edge) * n2 / (2.0 * area2 * edge_norm);
        u4 = -dot(ve1 - ve3, edge) * n1 / (2.0 * area1 * edge_norm) - dot(ve2 - ve3, edge) * n2 / (2.0 * area2 * edge_norm);
        
        forces[3 * vert_1 + 0] += C * u1.x;
        forces[3 * vert_1 + 1] += C * u1.y;
        forces[3 * vert_1 + 2] += C * u1.z;
        
        forces[3 * vert_2 + 0] += C * u2.x;
        forces[3 * vert_2 + 1] += C * u2.y;
        forces[3 * vert_2 + 2] += C * u2.z;
        
        forces[3 * vert_3 + 0] += C * u3.x;
        forces[3 * vert_3 + 1] += C * u3.y;
        forces[3 * vert_3 + 2] += C * u3.z;
        
        forces[3 * vert_4 + 0] += C * u4.x;
        forces[3 * vert_4 + 1] += C * u4.y;
        forces[3 * vert_4 + 2] += C * u4.z;
    }
}


double ForcesCalculator::calculate_theta(const Vector3D& ve1, 
                                         const Vector3D& ve2, 
                                         const Vector3D& ve3, 
                                         const Vector3D& ve4) const
{   
    Vector3D edge = ve4 - ve3;
    Vector3D e_norm = edge / edge.length();

    Vector3D N1 = cross(ve3 - ve1, ve4 - ve1);
    Vector3D n1 = N1 / N1.length();
    Vector3D N2 = cross(ve4 - ve2, ve3 - ve2);
    Vector3D n2 = N2 / N2.length();

    Vector3D cross_n1n2 = cross(n1, n2);
    int sign = SIGN( dot( cross_n1n2, e_norm) );

    double sin_theta = sign * cross_n1n2.length();
    
    return asin(sin_theta);
}


void ForcesCalculator::evaluate_pressure(const std::vector<double>& xyz, 
                          std::vector<double>& forces,
                          const std::vector<element>& elements,
                          const std::vector<object_map> vs_map,
                          const std::vector<double> turgors,
                          const int num_shells) const
{
    std::vector<Vector3D> cms;
    std::vector<int> vertex_counter;
    
    for (int i = 0; i < num_shells; i++)
    {
        cms.push_back(Vector3D(0,0,0));
        vertex_counter.push_back(0);
    }
    
    int cell_id;
    for (uint i = 0; i < xyz.size()/3; i++)
    {
        cell_id = vs_map[i].cell_id;
        cms[cell_id].x += xyz[3*i + 0];
        cms[cell_id].y += xyz[3*i + 1];
        cms[cell_id].z += xyz[3*i + 2];
        vertex_counter[cell_id]++;
    }
    
    for (uint i = 0; i < cms.size(); i++)
    {
        cms[i] /= vertex_counter[i];
    }
    
    
    double x1, y1, z1;
    double x2, y2, z2;
    double x3, y3, z3;
    double turgor;
    
    int vert_a, vert_b, vert_c;
    Vector3D va, vb, vc;
    Vector3D cm;
    element el;
    
    for (uint i = 0; i < elements.size(); i++)
    {
        el = elements[i];
        vert_a = el.ia;
        vert_b = el.ib;
        vert_c = el.ic;
        
        if (vs_map[vert_a].cell_id == vs_map[vert_b].cell_id && vs_map[vert_a].cell_id == vs_map[vert_c].cell_id)
        {
            cell_id = vs_map[i].cell_id;
            cm = cms[cell_id];
        }
        else
        {
            //PRINT ERROR AND ABORT
        }
        
        x1 = xyz[3 * vert_a + 0];
        y1 = xyz[3 * vert_a + 1];
        z1 = xyz[3 * vert_a + 2];
        
        x2 = xyz[3 * vert_b + 0];
        y2 = xyz[3 * vert_b + 1];
        z2 = xyz[3 * vert_b + 2];
        
        x3 = xyz[3 * vert_c + 0];
        y3 = xyz[3 * vert_c + 1];
        z3 = xyz[3 * vert_c + 2];
        
        va = Vector3D(x1, y1, z1);
        vb = Vector3D(x2, y2, z2);
        vc = Vector3D(x3, y3, z3);
        
        turgor = turgors[cell_id];
        
        Vector3D fa = turgor * calculate_dV(va, vb, vc, cm);
        Vector3D fb = turgor * calculate_dV(vb, vc, va, cm);
        Vector3D fc = turgor * calculate_dV(vc, va, vb, cm);
        
        
        forces[3 * vert_a + 0] += fa.x;
        forces[3 * vert_a + 1] += fa.y;
        forces[3 * vert_a + 2] += fa.z;
        
        forces[3 * vert_b + 0] += fb.x;
        forces[3 * vert_b + 1] += fb.y;
        forces[3 * vert_b + 2] += fb.z;
        
        forces[3 * vert_c + 0] += fc.x;
        forces[3 * vert_c + 1] += fc.y;
        forces[3 * vert_c + 2] += fc.z;  
    }
    
}

Vector3D ForcesCalculator::calculate_dV(const Vector3D& va,
                                        const Vector3D& vb,
                                        const Vector3D& vc,
                                        const Vector3D& vd) const
{
    Vector3D BD = vb - vd;
    Vector3D CD = vc - vd;
    Vector3D f = cross(BD, CD) / 6;

    return f * Tetrahedron::volumeSgn(va, vb, vc, vd);
}