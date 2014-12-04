#include "Aspherity.h"

Aspherity::Aspherity() {}

Aspherity::Aspherity(const Aspherity& orig) {}

Aspherity::~Aspherity() {}

double Aspherity::calcAspherity(Cell cell)
{
    double av_radius = 0.0;
    double sq_av_radius = 0.0;
    double sq_sum = 0.0;
    cell.calcCM();
    Vector3D cell_cm = cell.cm;
    
    for (int i = 0; i < cell.numberOfVerts(); i++)
    {
        av_radius += (cell.vertices[i].xyz - cell_cm).length();
    }
    av_radius /= cell.numberOfVerts();
    sq_av_radius = av_radius * av_radius;
    
    double res;
    for (int i = 0; i < cell.numberOfVerts(); i++)
    {
        res = (cell.vertices[i].xyz - cell_cm).length() - av_radius;
        sq_sum += res * res;
    }
    
    sq_sum /= sq_av_radius;
    sq_sum /= cell.numberOfVerts();
    
    return sq_sum;
}

