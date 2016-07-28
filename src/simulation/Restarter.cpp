#include "Restarter.h"

utils::Logger Restarter::restarter_logs("restarter");

Restarter::Restarter(std::string tf, std::string lff) : topologyFile(tf), lastFrameFile(lff)
{}

Restarter::Restarter(const Restarter& orig) : topologyFile(orig.topologyFile), lastFrameFile(orig.lastFrameFile) {}

Restarter::~Restarter()
{
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
        os << "NUMCELLS " << cells.size() << ' ' << model_t << "\n";

        for (uint i = 0; i < cells.size(); i++)
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

void Restarter::readTopologyFile(std::vector<Cell>& cells) const
{
    std::pair<int,std::string>  nc_mtype = getNumberOfCells();

    for (int i = 0; i < nc_mtype.first; i++)
    {
        Cell newCell;
        cells.push_back(newCell);
    }
    
    
    std::cout << cells.size() << " " << nc_mtype.first << std::endl;
    
    for (int i = 0; i < nc_mtype.first; i++)
    {
        std::cout << nc_mtype.second << std::endl;
        if ( nc_mtype.second.compare("fem") == 0 )
        {
            cells[i].fem_flag = true;
            std::cout << "FEM TRUE" << std::endl;
        }
    }
    
    for (int i = 0; i < nc_mtype.first; i++)
    {
        initCell(cells, i);
    }
    
}

std::pair<int,std::string> Restarter::getNumberOfCells() const
{  
    std::ifstream os;
    os.open(topologyFile, std::ifstream::in);
    std::vector<std::string> list;
    std::string line;
    

    std::pair<int,std::string> line_pair;//(-1,"NULL");
    
    if ( os.is_open() )
    {

        while ( std::getline (os, line) )
        {
            if ( line.find("NUMCELLS") == 0 )
            {
                std::vector<std::string> pairs = split(line, ' ');
                
                int n_cells = std::stoi(pairs[ 1 ].c_str(), NULL);
                std::string model_type(pairs[2]);
                
                //std::cout << typeid(model_type).name()<< std::endl;
                
                line_pair = std::pair<int,std::string>(n_cells, model_type);
            }
        }
    }
    else
    {
       // print error  
    }

    os.close();
    return line_pair;
}

void Restarter::initCell(std::vector<Cell>& cells, int cix) const
{  
    std::ifstream os;
    os.open(topologyFile, std::ifstream::in);
    std::vector<std::string> list;
    std::string line;
    
    if ( os.is_open() )
    {

        while ( std::getline (os, line) )
        {
            if ( line.find("CELL ") == 0 )
            {
                std::vector<std::string> pairs = split(line, ' ');
                
                int cell_id = std::stoi(pairs[ 1 ].c_str(), NULL);
                if (cell_id == cix)
                {
                    cells[cix].cell_id  = std::stoi(pairs[ 1 ].c_str(), NULL);
                    cells[cix].number_v  = std::stoi(pairs[ 2 ].c_str(), NULL);
                    cells[cix].number_t  = std::stoi(pairs[ 3 ].c_str(), NULL);
                    cells[cix].number_s  = std::stoi(pairs[ 4 ].c_str(), NULL);
                    
                    cells[cix].params.vertex_r = strtod(pairs[ 5 ].c_str(), NULL);
                    cells[cix].params.ecc = strtod(pairs[ 6 ].c_str(), NULL);
                    cells[cix].params.nu = strtod(pairs[ 7 ].c_str(), NULL);
                    cells[cix].params.dp = strtod(pairs[ 8 ].c_str(), NULL);
                    cells[cix].params.init_r = strtod(pairs[ 9 ].c_str(), NULL);
                    cells[cix].params.vol_c = strtod(pairs[ 10 ].c_str(), NULL);
                    cells[cix].nRT = strtod(pairs[ 11 ].c_str(), NULL);
                    cells[cix].V0 = strtod(pairs[ 12 ].c_str(), NULL);
                }
            }
        }
    }
    else
    {
       // print error  
    }

    os.close();
}


void Restarter::addVertices(std::vector<Cell>& cells, int cix) const
{
    std::ifstream os;
    os.open(topologyFile, std::ifstream::in);
    std::vector<std::string> list;
    std::string line;
    
    if ( os.is_open() )
    {
        while ( std::getline (os, line) )
        {
            if ( line.find("CELLVERTEX ") == 0 )
            {
                std::vector<std::string> pairs = split(line, ' ');
                
                int cell_id = std::stoi(pairs[ 1 ].c_str(), NULL);
                if (cell_id == cix)
                {
                    int v_id = std::stoi(pairs[ 2 ].c_str(), NULL);
                    //cells[cix].vertices[v_id];
                }
            }
        }
    }
    else
    {
       // print error  
    }

    os.close();    
}
