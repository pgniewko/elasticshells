#include "Element.h"

double cot(double x)
{
    return 1.0 / tan(x);

}

Element::Element() : Element(-1, -1, -1)
{
}

Element::Element(int a, int b, int c) : ia(a), ib(b), ic(c), my_id(-1)
{
    for (int i = 0; i < 3; i++)
    {
        an[i] = 0.0;
        L2[i] = 0.0;
        ki[i] = 0.0;
        ci[i] = 0.0;
    }
}

Element::Element(const Element& orig) : ia(orig.ia),
    ib(orig.ib),
    ic(orig.ic),
    my_id(orig.my_id)
{
    for (int i = 0; i < 3; i++)
    {
        an[i] = orig.an[i];
        L2[i] = orig.L2[i];
        ki[i] = orig.ki[i];
        ci[i] = orig.ci[i];
    }
}

Element::~Element() {}

void Element::set_id(int idx)
{
    my_id = idx;
}

int Element::get_id() const
{
    return my_id;
}

double Element::area(const std::vector<Vertex>& vs) const
{
    if (ia != -1 && ib != -1 && ic != -1)
    {
        Triangle t(vs[ia].r_c, vs[ib].r_c, vs[ic].r_c);
        return t.area();
    }
    else
    {
        return 0.0;
    }
}

double Element::area(const std::vector<Vertex>& vs, const Vector3D cm, double eps) const
{
    if (eps == 0)
    {
        return area(vs);
    }

    if (ia != -1 && ib != -1 && ic != -1)
    {
        Vector3D ta = vs[ia].r_c - cm;
        Vector3D tb = vs[ib].r_c - cm;
        Vector3D tc = vs[ic].r_c - cm;
        double n_ta = ta.length() + eps;
        double n_tb = tb.length() + eps;
        double n_tc = tc.length() + eps;
        ta.set_length(n_ta);
        tb.set_length(n_tb);
        tc.set_length(n_tc);
        Triangle t(ta, tb, tc);
        return t.area();
    }
    else
    {
        return 0.0;
    }
}

Vector3D Element::normal(const std::vector<Vertex>& vs) const
{
    Vector3D ta = vs[ia].r_c;
    Vector3D tb = vs[ib].r_c;
    Vector3D tc = vs[ic].r_c;
    Triangle t(ta, tb, tc);
    Vector3D normal = t.normal();
    return normal;
}

void Element::set_l2(const std::vector<Vertex>& vs)
{
    L2[0] = (vs[ib].r_c - vs[ic].r_c).length_sq();
    L2[1] = (vs[ia].r_c - vs[ic].r_c).length_sq();
    L2[2] = (vs[ia].r_c - vs[ib].r_c).length_sq();
}

void Element::set_an(const std::vector<Vertex>& vs)
{
    // MAKE SURE THAT the angle is between 0-180
    Vector3D ca = vs[ic].r_c - vs[ia].r_c;
    Vector3D ba = vs[ib].r_c - vs[ia].r_c;
    an[0] = ca.angle(ba);

    Vector3D ab = vs[ia].r_c - vs[ib].r_c;
    Vector3D cb = vs[ic].r_c - vs[ib].r_c;
    an[1] = ab.angle(cb);

    Vector3D ac = vs[ia].r_c - vs[ic].r_c;
    Vector3D bc = vs[ib].r_c - vs[ic].r_c;
    an[2] = ac.angle(bc);
}

void Element::set_ki(const std::vector<Vertex>& vs, const double& E, const double& nu, const double& t)
{
    double Ap = area(vs);
    ki[0] = E * t * (2.0 * cot(an[0]) * cot(an[0]) + 1.0 - nu) / (16.0 * Ap * (1.0 - nu * nu));
    ki[1] = E * t * (2.0 * cot(an[1]) * cot(an[1]) + 1.0 - nu) / (16.0 * Ap * (1.0 - nu * nu));
    ki[2] = E * t * (2.0 * cot(an[2]) * cot(an[2]) + 1.0 - nu) / (16.0 * Ap * (1.0 - nu * nu));

}

void Element::set_ci(const std::vector<Vertex>& vs, const double& E, const double& nu, const double& t)
{
    double Ap = area(vs);
    ci[0] = E * t * (2.0 * cot(an[1]) * cot(an[2]) + nu - 1.0 ) / (16.0 * Ap * (1.0 - nu * nu));
    ci[1] = E * t * (2.0 * cot(an[0]) * cot(an[2]) + nu - 1.0 ) / (16.0 * Ap * (1.0 - nu * nu));
    ci[2] = E * t * (2.0 * cot(an[0]) * cot(an[1]) + nu - 1.0 ) / (16.0 * Ap * (1.0 - nu * nu));
}

void Element::set_params(const std::vector<Vertex>& vs, const double E, const double nu, const double t)
{
    set_l2(vs);
    set_an(vs);
    set_ki(vs, E, nu, t);
    set_ci(vs, E, nu, t);
}

std::ostream& operator<< (std::ostream& out, const Element& vt)
{
    out << vt.my_id << ' ' << vt.ia << ' ' << vt.ib << ' ' << vt.ic << ' ';
    out << vt.an[0] << ' ' << vt.an[1] << ' '  << vt.an[2] << ' ';
    out << vt.L2[0] << ' ' << vt.L2[1] << ' '  << vt.L2[2] << ' ';
    out << vt.ki[0] << ' ' << vt.ki[1] << ' '  << vt.ki[2] << ' ';
    out << vt.ci[0] << ' ' << vt.ci[1] << ' '  << vt.ci[2] << ' ';
    return out;
}