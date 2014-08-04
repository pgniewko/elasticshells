#include "Triangulation.h"

Triangulation::Triangulation() {}

Triangulation::Triangulation(const Triangulation& orig) {}

Triangulation::~Triangulation() {}

void Triangulation::saveTriangulatedSurface(const char* filename, bool wrap)
{
    ofstream os(filename);
    
    int index = 1;
    os << 3*(int)tris.size() << "\n" ;
    for(list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i) // go through all triangles
    {
        if (wrap)
        {
//      write
        }
        else
        {
           os << "H" << index++ << " "<< i->a.x << " " << i->a.y << " " << i->a.z << "\n";
           os << "H" << index++ << " "<< i->b.x << " " << i->b.y << " " << i->b.z << "\n";
           os << "H" << index++ << " "<< i->c.x << " " << i->c.y << " " << i->c.z << "\n";
        }
    }

    os.close();
}

void Triangulation::saveRenderingScript(const char* filename, const char* cellsfile)
{
    ofstream os(filename);
    
    os << "from pymol.cgo import *\n";
    os << "from pymol import cmd \n\n";
    os << "cmd.do(\"load " << cellsfile << ", cells\")\n";
    os << "cmd.do(\"hide all\")\n";
    os << "cmd.do(\"set sphere_color, tv_red\")\n";
    os << "cmd.do(\"set line_color, marine\")\n";
    os << "cmd.do(\"show spheres\")\n";
    os << "cmd.do(\"alter elem h, vdw=0.1\")\n";
    os << "cmd.do(\"rebuild\")\n";


    int index = 1;
    for(std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i)
    {
        os << "cmd.do(\"bond /cells///UNK`/H"<< index+0 << ", /cells///UNK`/H" << index+1 << "\")\n";
        os << "cmd.do(\"bond /cells///UNK`/H"<< index+1 << ", /cells///UNK`/H" << index+2 << "\")\n";
        os << "cmd.do(\"bond /cells///UNK`/H"<< index+2 << ", /cells///UNK`/H" << index+0 << "\")\n";
        index = index + 3;
    }
    
    os << "cmd.do(\"show lines\")\n";
    os << "cmd.do(\"bg white\")\n";

    os.close();
}