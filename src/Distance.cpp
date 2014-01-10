#include "Distance.h"

Distance::Distance() 
{
    lx = 0.0;
    ly = 0.0;
    lz = 0.0;
}

Distance::Distance(double l_x, double l_y, double l_z, Vector3D o) 
{
    this->set(l_x, l_y, l_z, o);    
}

Distance::~Distance() 
{
}

void Distance::set(double l_x, double l_y, double l_z, Vector3D o)
{
    lx = l_x;
    ly = l_y;  
    lz = l_z;  
    ilx = (lx? 1.0 / lx : 0.0);
    ily = (ly? 1.0 / ly : 0.0);
    ilz = (lz? 1.0 / lz : 0.0);
    O = o;
}

Vector3D Distance::delta(const Vector3D& p1) const
{
    Vector3D delta = p1 - O;
    delta.x -= lx * rint( ilx * delta.x );
    delta.y -= ly * rint( ily * delta.y );
    delta.z -= lz * rint( ilz * delta.z );
    return delta;
}

Vector3D Distance::delta(const Vector3D& p1, const  Vector3D& p2) const
{
    Vector3D delta = p1 - p2;
    delta.x -= lx * rint( ilx * delta.x );
    delta.y -= ly * rint( ily * delta.y );
    delta.z -= lz * rint( ilz * delta.z );
    return delta;
}
	
Vector3D Distance::image(const Vector3D& p1) const 
{
    Vector3D delta = p1 - O;
    delta.x -= lx * rint( ilx * delta.x );
    delta.y -= ly * rint( ily * delta.y );
    delta.z -= lz * rint( ilz * delta.z );
    return delta + O;
}