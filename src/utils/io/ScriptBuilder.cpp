#include "ScriptBuilder.h"

utils::Logger ScriptBuilder::scriptbuilder_logs("scriptbuilder");

ScriptBuilder::ScriptBuilder(std::string rs, std::string ss, std::string tf, std::string sx) :
    script(rs), surfaceScript(ss), trajfile(tf), stress_script(sx) {}

ScriptBuilder::ScriptBuilder(const ScriptBuilder& orig) :
    script(orig.script), surfaceScript(orig.surfaceScript), trajfile(orig.trajfile), stress_script(orig.stress_script) {}

ScriptBuilder::~ScriptBuilder() {}

void ScriptBuilder::saveRenderScript(const std::vector<Cell>& cells, const Box& box, bool boxFlag, double rv) const
{
    std::ofstream os;
    os.open(script);

    if ( os.is_open() )
    {
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
            printBox(box, os);
        }

        int lastCellIndex = 0;

        for (uint i = 0; i < cells.size(); i++)
        {
            for (int j = 0; j < cells[i].getNumberVertices(); j++)
            {
                std::string strindex = new_base_index( lastCellIndex +  cells[i].vertices[j].getId());
                os << "cmd.do(\"select " << strindex << ", name " << strindex << "\")\n";
            }

            lastCellIndex += cells[i].getNumberVertices();
        }

        lastCellIndex = 0;

        for (uint i = 0; i < cells.size(); i++)
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
    else
    {
        scriptbuilder_logs << utils::LogLevel::WARNING << "Can not open file:" <<  script << "\n";
    }
}

void ScriptBuilder::printBox(const Box& box, std::ofstream& os) const
{
    if ( os.is_open() )
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
    else
    {
        scriptbuilder_logs << utils::LogLevel::WARNING << "No open file for box script.\n";
    }
}