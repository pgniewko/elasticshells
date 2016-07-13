#include "Restarter.h"

utils::Logger Restarter::restarter_logs("restarter");

Restarter::Restarter(std::string tf, std::string lff) :topologyFile(tf), lastFrameFile(lff)
{}

Restarter::Restarter(const Restarter& orig) : topologyFile(orig.topologyFile), lastFrameFile(orig.lastFrameFile) {}

Restarter::~Restarter() {
}


void Restarter::saveTopologyFile(const std::vector<Cell>& cells) const
{
    std::ofstream os;
    os.open(topologyFile);
    if ( os.is_open() )
    {
        os << "NUMCELLS " << cells.size() << "\n";
        for (uint i=0; i<cells.size();i++)
        {
            os << cells[i] << "\n";
            
        }
        
        os.close();
    }
    else
    {
        restarter_logs << utils::LogLevel::WARNING << "Could not open file:" << topologyFile << "\n";
    }
}
