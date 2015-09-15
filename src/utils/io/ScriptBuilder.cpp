#include "ScriptBuilder.h"

extern const char names[];

ScriptBuilder::ScriptBuilder(char* rs, char* ss, char* tf, char* sx)
{
    script = rs;
    surfaceScript = ss;
    trajfile = tf;
    stress_script = sx;
}

ScriptBuilder::ScriptBuilder(const ScriptBuilder& orig) {}

ScriptBuilder::~ScriptBuilder() {}

void ScriptBuilder::setDrawBox(bool db)
{
    drawBox = db;
}


void ScriptBuilder::saveSurfaceScript(const std::vector<Cell>& cells)
{
    std::ofstream os;
    os.open(surfaceScript);
    os << "from pymol.cgo import *\n";
    os << "from pymol import cmd \n\n";
    os << "def compute_normal(x1, y1, z1, x2, y2, z2, x3, y3, z3):\n\n";
    os << "  nx = (y2-y1)*(z3-z2) - (z2-z1)*(y3-y2)\n";
    os << "  ny = (z2-z1)*(x3-x2) - (x2-x1)*(z3-z2)\n";
    os << "  nz = (x2-x1)*(y3-y2) - (y2-y1)*(x3-x2)\n\n";
    os << "  return (nx,ny,nz)\n\n\n";
    os << "def draw_plane_cgo(name, apex1, apex2, apex3, color=(1,1,1)):\n";
    os << "  x1,y1,z1 = map(float,apex1)\n";
    os << "  x2,y2,z2 = map(float,apex2)\n";
    os << "  x3,y3,z3 = map(float,apex3)\n";
    os << "  if type(color) == type(''):\n";
    os << "    color = map(float,color.replace('(','').replace(')','').split(','))\n\n";
    os << "  # Compute the normal vector for the triangle\n";
    os << "  normal1 = compute_normal(x1, y1, z1, x2, y2, z2, x3, y3, z3)\n\n";
    os << "  # Create the CGO objects\n";
    os << "  obj = [\n";
    os << "    BEGIN, TRIANGLE_STRIP,\n\n";
    os << "    COLOR, color[0], color[1], color[2],\n";
    os << "    NORMAL, normal1[0], normal1[1], normal1[2],\n";
    os << "    NORMAL, -normal1[0], -normal1[1], -normal1[2],\n";
    os << "    VERTEX, x1, y1, z1,\n";
    os << "    VERTEX, x2, y2, z2,\n";
    os << "    VERTEX, x3, y3, z3,\n\n";
    os << "    END\n";
    os << "  ]\n\n";
    os << "  # Display them\n";
    os << "  cmd.load_cgo(obj, name)\n\n";
    os << "def draw_plane(name,atom1='(pk1)',atom2='(pk2)',atom3='(pk3)',color=(1,1,1)):\n";
    os << "# get coordinates from atom selections\n";
    os << "  coor1 = cmd.get_model(atom1).atom[0].coord\n";
    os << "  coor2 = cmd.get_model(atom2).atom[0].coord\n";
    os << "  coor3 = cmd.get_model(atom3).atom[0].coord\n";
    os << "  draw_plane_cgo(name,coor1,coor2,coor3,color)\n\n";
    os << "cmd.extend(\"draw_plane\", draw_plane)\n\n\n";
    int lastCellIndex = 0;
    int faceCounter = 0;
    int idxa, idxb, idxc;

    for (unsigned int i = 0; i < cells.size(); i++)
    {
        for (int j = 0; j < cells[i].getNumberTriangles(); j++)
        {
            idxa = lastCellIndex + cells[i].triangles[j].ia;
            idxb = lastCellIndex + cells[i].triangles[j].ib;
            idxc = lastCellIndex + cells[i].triangles[j].ic;
            std::string strindex1 = new_base_index( idxa );
            std::string strindex2 = new_base_index( idxb );
            std::string strindex3 = new_base_index( idxc );

            os << "cmd.do(\"draw_plane \\\"face" << faceCounter << "\\\", " << strindex1 << ", " <<  strindex2 << ", " << strindex3 << ", (0.8, 0.8, 0.8) \")\n";
            faceCounter++;
        }

        lastCellIndex += cells[i].getNumberVertices();
    }

    os.close();
}


void ScriptBuilder::saveRenderScript(const std::vector<Cell>& cells, const Box& box, bool boxFlag, double rv)
{
    std::ofstream os;
    os.open(script);
    os << "from pymol.cgo import *\n";
    os << "from pymol import cmd \n\n";
    os << "cmd.do(\"load " << trajfile << ", cells\")\n";
    os << "cmd.do(\"select rawdata, all\")\n";
    os << "cmd.do(\"unbond rawdata, rawdata\")\n";
    os << "cmd.do(\"hide all\")\n";
    os << "cmd.do(\"set sphere_color, tv_red\")\n";
    os << "cmd.do(\"set line_color, marine\")\n";
    os << "cmd.do(\"show spheres\")\n";
    os << "cmd.do(\"alter cells, vdw=" << rv << "\")\n";
    os << "cmd.do(\"rebuild\")\n\n";

    if (boxFlag)
    {
        printBox(os, box);
    }

    int lastCellIndex = 0;

    for (unsigned int i = 0; i < cells.size(); i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            std::string strindex = new_base_index( lastCellIndex +  cells[i].vertices[j].getId());
            os << "cmd.do(\"select " << strindex << ", name " << strindex << "\")\n";
        }

        lastCellIndex += cells[i].getNumberVertices();
    }

    lastCellIndex = 0;

    for (unsigned int i = 0; i < cells.size(); i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            for (int k = 0; k < cells[i].vertices[j].numBonded; k++)
            {
                std::string strindex1 = new_base_index( lastCellIndex +  cells[i].vertices[j].getId());
                std::string strindex2 = new_base_index( lastCellIndex +  cells[i].vertices[j].bondedVerts[k]);
                os << "cmd.do(\"bond " << strindex1 << ", " << strindex2 << "\")\n";
            }
        }

        lastCellIndex += cells[i].getNumberVertices();
    }

    os << "cmd.do(\"show lines\")\n";
    os << "cmd.do(\"bg white\")\n\n";
    os.close();
}

void ScriptBuilder::printBox(std::ofstream& os, const Box& box)
{
    os << "class Box(object):\n";
    os << "  def __init__ (self, x, y, z, linewidth, color):\n";
    os << "    lw = linewidth\n";
    os << "    c1 = color[0]\n";
    os << "    c2 = color[1]\n";
    os << "    c3 = color[2]\n";
    os << "    self.box = [\n";
    os << "    LINEWIDTH, float(lw), BEGIN, LINES,\n";
    os << "    COLOR,  color[0], color[1], color[2],\n";
    os << "    VERTEX, x[0], y[0], z[0],\n";
    os << "    VERTEX, x[1], y[0], z[0],\n";
    os << "    VERTEX, x[0], y[0], z[0],\n";
    os << "    VERTEX, x[0], y[1], z[0],\n";
    os << "    VERTEX, x[0], y[0], z[0],\n";
    os << "    VERTEX, x[0], y[0], z[1],\n";
    os << "    VERTEX, x[1], y[1], z[1],\n";
    os << "    VERTEX, x[1], y[1], z[0],\n";
    os << "    VERTEX, x[1], y[1], z[1],\n";
    os << "    VERTEX, x[0], y[1], z[1],\n";
    os << "    VERTEX, x[1], y[1], z[1],\n";
    os << "    VERTEX, x[1], y[0], z[1],\n";
    os << "    VERTEX, x[0], y[0], z[1],\n";
    os << "    VERTEX, x[1], y[0], z[1],\n";
    os << "    VERTEX, x[0], y[1], z[0],\n";
    os << "    VERTEX, x[0], y[1], z[1],\n";
    os << "    VERTEX, x[0], y[1], z[0],\n";
    os << "    VERTEX, x[1], y[1], z[0],\n";
    os << "    VERTEX, x[0], y[0], z[1],\n";
    os << "    VERTEX, x[0], y[1], z[1],\n";
    os << "    VERTEX, x[1], y[0], z[0],\n";
    os << "    VERTEX, x[1], y[1], z[0],\n";
    os << "    VERTEX, x[1], y[0], z[0],\n";
    os << "    VERTEX, x[1], y[0], z[1],\n";
    os << "    END\n";
    os << "    ]\n\n\n";
    os << "B = Box(";
    os << "(" << -box.getX() << "," << box.getX() << "),";
    os << "(" << -box.getY() << "," << box.getY() << "),";
    os << "(" << -box.getZ() << "," << box.getZ() << "),";
    os << " 2.5, color=(0.0, 0.0, 0.0) )\n";
    os << "obj = B.box\n";
    os << "cmd.load_cgo(obj, \"box\", 1)\n\n\n";
}

void ScriptBuilder::saveStressScript(std::vector<Cell>& cells, const Box& box)
{
    std::ofstream os;
    os.open(stress_script);
    os << "from pymol.cgo import *\n";
    os << "from pymol import cmd \n\n";
    os << "class Line(object):\n";
    os << " def __init__ (self, x, y, z, linewidth, color):\n";
    os << "  lw = linewidth\n";
    os << "  c1 = color[0]\n";
    os << "  c2 = color[1]\n";
    os << "  c3 = color[2]\n";
    os << "  self.line = [\n";
    os << "    CYLINDER,  x[0], y[0], z[0], x[1], y[1], z[1], linewidth, c1, c2,c3,c1,c2,c3, \n";
    os << "    END\n";
    os << "  ]\n\n";
    os << "cmd.do(\"hide spheres\")\n";
    int counter = 1;
    double maxforce = 0.0;

    for (unsigned int i = 0; i < cells.size(); i++)
    {
        for (unsigned int j = 0; j < cells.size(); j++)
        {
            double forceij = 0.0;
            maxforce = std::max(maxforce, forceij);
        }
    }

    maxforce *= 2.0;

    for (unsigned int i = 0; i < cells.size(); i++)
    {
        for (unsigned int j = 0; j < cells.size(); j++)
        {
            double forceij = 0.0;
            //double forceij = cells[i].contactForceNew(cells[j], box);
            //cells[i].calcCM();
            Vector3D cmi = cells[i].getCm();
            //cells[j].calcCM();
            Vector3D cmj = cells[j].getCm();

            if (forceij > 0 and i > j)
            {
                double xi = cmi.x * box.getXstart() / box.getX();
                double xj = cmj.x * box.getXstart() / box.getX();
                double yi = cmi.y * box.getYstart() / box.getY();
                double yj = cmj.y * box.getYstart() / box.getY();
                double zi = cmi.z * box.getZstart() / box.getZ();
                double zj = cmj.z * box.getZstart() / box.getZ();
                os << "L = Line(";
                os << "(" << xi << "," << xj << "),";
                os << "(" << yi << "," << yj << "),";
                os << "(" << zi << "," << zj << "), ";
                //os << "(" << forceij / maxforce << "), ";
                os << "(" << forceij << "), ";
                os << "color=(1.0, 0.0, 0.0) )\n";
                os << "obj = L.line\n";
                os << "cmd.load_cgo(obj, \"line" << counter << "\", 1)\n";
                counter++;
            }
        }
    }

    os.close();
}