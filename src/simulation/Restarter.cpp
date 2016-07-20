#include "Restarter.h"

utils::Logger Restarter::restarter_logs("restarter");

Restarter::Restarter(std::string tf, std::string lff) :topologyFile(tf), lastFrameFile(lff)
{}

Restarter::Restarter(const Restarter& orig) : topologyFile(orig.topologyFile), lastFrameFile(orig.lastFrameFile) {}

Restarter::~Restarter() {
}

int Restarter::getTotalVertices(const std::vector<Cell>& cells) const
{
    int totalnumber = 0;

    for (uint i = 0; i < cells.size(); i++)
    {
        totalnumber += cells[i].getNumberVertices();
    }

    return totalnumber;
}

void Restarter::saveTopologyFile(const std::vector<Cell>& cells, char* model_t) const
{
    std::ofstream os;
    os.open(topologyFile);
    if ( os.is_open() )
    {
        os << "NUMCELLS " << cells.size() << ' ' << model_t<< "\n";
        for (uint i = 0; i< cells.size(); i++)
        {
            os << cells[i];
            
        }
        
        int lastCellIndex = 0;
        
        for (uint i = 0; i < cells.size(); i++)
        {
            for (int j = 0; j < cells[i].getNumberVertices(); j++)
            {
                os << "VMAP " <<  new_base_index( lastCellIndex +  cells[i].vertices[j].getId() ) << ' ' << i << ' ' << j << "\n";
            }

            lastCellIndex += cells[i].getNumberVertices();
        }
        
        os.close();
    }
    else
    {
        restarter_logs << utils::LogLevel::WARNING << "Could not open file:" << topologyFile << "\n";
    }
}

void Restarter::saveLastFrame(const std::vector<Cell>& cells) const
{
    XyzTraj lf_xyz(lastFrameFile, "NULL");
    lf_xyz.open_traj();
    lf_xyz.save_traj(cells, getTotalVertices(cells));
    lf_xyz.close_traj();
}
