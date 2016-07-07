#include "Triangulation.h"

Triangulation::Triangulation() {}

Triangulation::Triangulation(const Triangulation& orig) {}

Triangulation::~Triangulation() {}

void Triangulation::saveTriangulatedSurface(const char* filename, bool wrap)
{
    std::ofstream os(filename);
    int index = 1;
    os << 3 * (int)tris.size() << "\n" ;

    for (std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i) // go through all triangles
    {
        if (wrap)
        {
//      write
        }
        else
        {
            std::string strindex = new_base_index( index );
            os << strindex << " " << i->a.x << " " << i->a.y << " " << i->a.z << "\n";
            index++;
            strindex = new_base_index( index );
            os << strindex << " " << i->b.x << " " << i->b.y << " " << i->b.z << "\n";
            index++;
            strindex = new_base_index( index );
            os << strindex << " " << i->c.x << " " << i->c.y << " " << i->c.z << "\n";
        }
    }

    os.close();
}

void Triangulation::saveRenderingScript(const char* filename, const char* cellsfile)
{
    std::ofstream os(filename);
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

    for (std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i)
    {
        std::string strindex1 = new_base_index( index + 0 );
        std::string strindex2 = new_base_index( index + 1 );
        std::string strindex3 = new_base_index( index + 2 );
        os << "cmd.do(\"bond /cells///UNK`/" << strindex1 << ", /cells///UNK`/" << strindex2 << "\")\n";
        os << "cmd.do(\"bond /cells///UNK`/" << strindex1 << ", /cells///UNK`/" << strindex3 << "\")\n";
        os << "cmd.do(\"bond /cells///UNK`/" << strindex2 << ", /cells///UNK`/" << strindex3 << "\")\n";
        index += 3;
    }

    os << "cmd.do(\"show lines\")\n";
    os << "cmd.do(\"bg white\")\n";
    os.close();
}