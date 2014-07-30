#include "Cell.h"

Cell::Cell(double x, double y, double z)
{
    r = Vector3D(x, y, z);
}

Vector3D Cell::reset_cforce()
{
    this->cf = Vector3D();
    return this->cf;
}

//Cell::Cell(const Cell& orig)
//{
//}

//Cell::~Cell()
//{
//}
