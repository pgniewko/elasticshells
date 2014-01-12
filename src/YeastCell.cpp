#include "YeastCell.h"

using namespace std;

YeastCell::YeastCell()
{
    params = new (nothrow) double[1];

    if (params == 0)
    {
        cout << "Error: memory could not be allocated";
    }
}

YeastCell::YeastCell(double r)
{
    params = new (nothrow) double[1];

    if (params == 0)
    {
        cout << "Error: memory could not be allocated";
    }

    params[0] = r;
}

YeastCell::YeastCell(const YeastCell& orig)
{
}

YeastCell::~YeastCell()
{
    delete[] params;
}

double YeastCell::calc_volume()
{
    double v = (4.0 * M_PI / 3.0) * params[0] * params[0] * params[0];
    v = params[0] * params[0] * params[0];
    return (v);
}

//void YeastCell::set_coor(double x, double y , double z)
//{
//
//}

void YeastCell::set_params(double* arr, int i)
{
//    cout << "uhu" << arr[0] << endl;
    params[0] = arr[0];
}

