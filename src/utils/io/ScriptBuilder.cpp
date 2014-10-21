#include "ScriptBuilder.h"

ScriptBuilder::ScriptBuilder(char* rs, char* ss, char* tf) : names( {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'})
{
    script = rs;
    surfaceScript = ss;
    trajfile = tf;
}

ScriptBuilder::ScriptBuilder(const ScriptBuilder& orig) {}

ScriptBuilder::~ScriptBuilder() {}

void ScriptBuilder::setDrawBox(bool db)
{
    drawBox = db;
}


void ScriptBuilder::saveSurfaceScript(vector<Cell>& cells)
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
    int name1Ix;
    int atom1Ix;
    int name2Ix;
    int atom2Ix;
    int name3Ix;
    int atom3Ix;
    char cA;
    char cB;
    char cC;

    for (unsigned int i = 0; i < cells.size(); i++)
    {
        for (int j = 0; j < cells[i].numberofFaces(); j++)
        {
            idxa = cells[i].triangles[j].ia + 1 + lastCellIndex;
            idxb = cells[i].triangles[j].ib + 1 + lastCellIndex;
            idxc = cells[i].triangles[j].ic + 1 + lastCellIndex;
            name1Ix = (int) idxa / 1000;
            atom1Ix = idxa % 1000;
            name2Ix = (int) idxb / 1000;
            atom2Ix = idxb % 1000;
            name3Ix = (int) idxc / 1000;;
            atom3Ix = idxc % 1000;
            cA = names[name1Ix];
            cB = names[name2Ix];
            cC = names[name3Ix];
            os << "cmd.do(\"draw_plane \\\"face" << faceCounter << "\\\", " << cA << "_" << atom1Ix << ", " << cB << "_" << atom2Ix << ", " << cC << "_" << atom3Ix << ", (0.8, 0.8, 0.8) \")\n";
            faceCounter++;
        }

        lastCellIndex += cells[i].numberofVertices();
    }

    os.close();
}


void ScriptBuilder::saveRenderScript(vector<Cell>& cells, Box& box, bool boxFlag)
{
    ofstream os;
    os.open(script);
    os << "from pymol.cgo import *\n";
    os << "from pymol import cmd \n\n";
    //if (boxFlag)
    //{
    //    printBox(os, box);
    //}
    os << "cmd.do(\"load " << trajfile << ", cells\")\n";
    os << "cmd.do(\"select rawdata, all\")\n";
    os << "cmd.do(\"unbond rawdata, rawdata\")\n";
    os << "cmd.do(\"hide all\")\n";
    os << "cmd.do(\"set sphere_color, tv_red\")\n";
    os << "cmd.do(\"set line_color, marine\")\n";
    os << "cmd.do(\"show spheres\")\n";
    os << "cmd.do(\"alter elem a, vdw=0.1\")\n";
    os << "cmd.do(\"alter elem b, vdw=0.1\")\n";
    os << "cmd.do(\"alter elem c, vdw=0.1\")\n";
    os << "cmd.do(\"alter elem d, vdw=0.1\")\n";
    os << "cmd.do(\"alter elem e, vdw=0.1\")\n";
    os << "cmd.do(\"alter elem f, vdw=0.1\")\n";
    os << "cmd.do(\"alter elem g, vdw=0.1\")\n";
    os << "cmd.do(\"alter elem h, vdw=0.1\")\n";
    os << "cmd.do(\"alter elem i, vdw=0.1\")\n";
    os << "cmd.do(\"alter elem j, vdw=0.1\")\n";
    os << "cmd.do(\"rebuild\")\n\n";

    if (boxFlag)
    {
        printBox(os, box);
    }

    char namesx[10] = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'};
    int name1Ix;
    int atom1Ix;
    int name2Ix;
    int atom2Ix;
    int iidx, jidx;
    int lastCellIndex = 0;

    for (unsigned int i = 0; i < cells.size(); i++)
    {
        for (int j = 0; j < cells[i].numberofVertices(); j++)
        {
            iidx = cells[i].vertices[j].getId() + 1 + lastCellIndex;
            name1Ix = (int) iidx / 1000;
            atom1Ix = iidx % 1000;
            os << "cmd.do(\"select " << namesx[name1Ix] << " " << atom1Ix << ", name " << namesx[name1Ix] << atom1Ix << "\")\n";
        }

        lastCellIndex += cells[i].numberofVertices();
    }

    lastCellIndex = 0;

    for (unsigned int i = 0; i < cells.size(); i++)
    {
        for (int j = 0; j < cells[i].numberofVertices(); j++)
        {
            iidx = cells[i].vertices[j].getId() + 1 + lastCellIndex;

            for (int k = 0; k < cells[i].vertices[j].nneigh; k++)
            {
                jidx = cells[i].vertices[j].neighbors[k] + 1 + lastCellIndex;
                name1Ix = (int) iidx / 1000;
                atom1Ix = iidx % 1000;
                name2Ix = (int) jidx / 1000;
                atom2Ix = jidx % 1000;
                os << "cmd.do(\"bond " << namesx[name1Ix] << "_" << atom1Ix << ", " << namesx[name2Ix] << "_" << atom2Ix << "\")\n";
            }
        }

        lastCellIndex += cells[i].numberofVertices();
    }

    os << "cmd.do(\"show lines\")\n";
    os << "cmd.do(\"bg white\")\n\n";
//    if (boxFlag)
//    {
//        os << "B = Box(";
//        os << "("<< -box.getX() << "," << box.getX() <<"),";
//        os << "("<< -box.getY() << "," << box.getY() <<"),";
//        os << "("<< -box.getZ() << "," << box.getZ() <<"),";
//        os << " 2.5, color=(0.0, 0.0, 0.0) )\n";
//        os << "obj = B.box\n";
//        os << "cmd.load_cgo(obj, \"box\", 1)\n";
//    }
    os.close();
}


void ScriptBuilder::printBox(ofstream& os,  Box& box)
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