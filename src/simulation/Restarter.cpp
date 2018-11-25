#include "Restarter.h"

utils::Logger Restarter::restarter_logs("restarter");

Restarter::Restarter(std::string tf, std::string lff) : topologyFile(tf), lastFrameFile(lff)
{}

Restarter::Restarter(const Restarter& orig) : topologyFile(orig.topologyFile), lastFrameFile(orig.lastFrameFile) {}

Restarter::~Restarter()
{
}

int Restarter::getTotalVertices(const std::vector<Shell>& shells) const
{
    int totalnumber = 0;

    for (uint i = 0; i < shells.size(); i++)
    {
        totalnumber += shells[i].getNumberVertices();
    }

    return totalnumber;
}

void Restarter::saveTopologyFile(const std::vector<Shell>& shells) const
{
    std::ofstream os;
    os.open(topologyFile);

    if ( os.is_open() )
    {
        os << "NUMSHELLS " << shells.size() << ' ' << "fem" << ' ' << (Shell::no_bending ? "true" : "false") << "\n";

        for (uint i = 0; i < shells.size(); i++)
        {
            os << shells[i];

        }

        int last_shell_index = 0;

        for (uint i = 0; i < shells.size(); i++)
        {
            for (int j = 0; j < shells[i].getNumberVertices(); j++)
            {
                os << "VMAP " <<  new_base_index( last_shell_index +  shells[i].vertices[j].getId() ) << ' ' << i << ' ' << j << "\n";
            }

            last_shell_index += shells[i].getNumberVertices();
        }

        os.close();
    }
    else
    {
        restarter_logs << utils::LogLevel::WARNING << "Could not open file:" << topologyFile << "\n";
    }
}

void Restarter::saveLastFrame(const std::vector<Shell>& shells, const Box& box) const
{
    XyzTraj lf_xyz(lastFrameFile, "NULL");
    lf_xyz.open_lf();
    lf_xyz.save_traj(shells, getTotalVertices(shells), box);
    lf_xyz.close_traj();
}

void Restarter::readTopologyFile(std::vector<Shell>& shells) const
{
    std::pair<int, std::string>  nc_mtype = get_number_of_shells();

    for (int i = 0; i < nc_mtype.first; i++)
    {
        Shell new_shell;
        shells.push_back(new_shell);
    }

    for (int i = 0; i < nc_mtype.first; i++)
    {
        initShell(shells, i);
    }

    for (int i = 0; i < nc_mtype.first; i++)
    {
        addVertices(shells, i);
    }

    for (int i = 0; i < nc_mtype.first; i++)
    {
        addVTriangles(shells, i);
    }

    for (int i = 0; i < nc_mtype.first; i++)
    {
        addBHinges(shells, i);
    }
}

std::pair<int, std::string> Restarter::get_number_of_shells() const
{
    std::ifstream os;
    os.open(topologyFile, std::ifstream::in);
    std::string line;


    std::pair<int, std::string> line_pair; //(-1,"NULL");

    if ( os.is_open() )
    {

        while ( std::getline (os, line) )
        {
            if ( line.find("NUMSHELLS") == 0 )
            {
                std::vector<std::string> pairs = split(line, ' ');

                int n_shells = std::stoi(pairs[ 1 ].c_str(), NULL);
                std::string model_type(pairs[2]);

                line_pair = std::pair<int, std::string>(n_shells, model_type);
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

void Restarter::initShell(std::vector<Shell>& shells, int cix) const
{
    std::ifstream os;
    os.open(topologyFile, std::ifstream::in);
    std::string line;

    if ( os.is_open() )
    {

        while ( std::getline (os, line) )
        {
            if ( line.find("SHELL ") == 0 )
            {
                std::vector<std::string> pairs = split(line, ' ');

                int shell_id = std::stoi(pairs[ 1 ].c_str(), NULL);

                if (shell_id == cix)
                {
                    shells[cix].shell_id  = std::stoi(pairs[ 1 ].c_str(), NULL);
                    shells[cix].number_v  = std::stoi(pairs[ 2 ].c_str(), NULL);
                    shells[cix].number_t  = std::stoi(pairs[ 3 ].c_str(), NULL);
                    shells[cix].number_s  = std::stoi(pairs[ 4 ].c_str(), NULL);

                    shells[cix].params.vertex_r = strtod(pairs[ 5 ].c_str(), NULL);
                    shells[cix].params.ecc = strtod(pairs[ 6 ].c_str(), NULL);
                    shells[cix].params.nu = strtod(pairs[ 7 ].c_str(), NULL);
                    shells[cix].params.dp = strtod(pairs[ 8 ].c_str(), NULL);
                    shells[cix].params.init_r = strtod(pairs[ 9 ].c_str(), NULL);
                    shells[cix].params.vol_c = strtod(pairs[ 10 ].c_str(), NULL);
                    shells[cix].nRT = strtod(pairs[ 11 ].c_str(), NULL);
                    shells[cix].V0 = strtod(pairs[ 12 ].c_str(), NULL);
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


void Restarter::addVertices(std::vector<Shell>& shells, int cix) const
{
    std::ifstream os;
    os.open(topologyFile, std::ifstream::in);
    std::string line;

    if ( os.is_open() )
    {
        while ( std::getline (os, line) )
        {
            if ( line.find("SHELLVERTEX ") == 0 )
            {
                std::vector<std::string> pairs = split(line, ' ');

                int shell_id = std::stoi(pairs[ 1 ].c_str(), NULL);

                if (shell_id == cix)
                {
                    int v_id = std::stoi(pairs[ 2 ].c_str(), NULL);
                    shells[cix].vertices[v_id].setId(v_id);
                    shells[cix].vertices[v_id].set_shell_id(cix);
                    shells[cix].vertices[v_id].numBonded = strtod(pairs[ 3 ].c_str(), NULL);
                    shells[cix].vertices[v_id].numTris = strtod(pairs[ 4 ].c_str(), NULL);

                    int start_ix = 4;

                    for (int i = 0; i < shells[cix].vertices[v_id].numBonded; i++)
                    {
                        shells[cix].vertices[v_id].bondedVerts[i] = std::stoi(pairs[ start_ix + 1 ].c_str(), NULL);
//                        shells[cix].vertices[v_id].r0[i]          = strtod(pairs[ start_ix + 2 ].c_str(), NULL);
//                        shells[cix].vertices[v_id].k0[i]          = strtod(pairs[ start_ix + 3 ].c_str(), NULL);
//                        start_ix += 3;
                        start_ix++;
                    }

                    start_ix++;

                    for (int i = 0; i < shells[cix].vertices[v_id].numTris; i++)
                    {
                        shells[cix].vertices[v_id].bondedTris[i] = std::stoi(pairs[ start_ix ].c_str(), NULL);
                        start_ix++;
                    }
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

void Restarter::addVTriangles(std::vector<Shell>& shells, int cix) const
{
    std::ifstream os;
    os.open(topologyFile, std::ifstream::in);
    std::string line;

    if ( os.is_open() )
    {
        while ( std::getline (os, line) )
        {
            if ( line.find("SHELLTRIANG ") == 0 )
            {
                std::vector<std::string> pairs = split(line, ' ');

                int shell_id = std::stoi(pairs[ 1 ].c_str(), NULL);

                if (shell_id == cix)
                {
                    int t_id = std::stoi(pairs[ 2 ].c_str(), NULL);
                    shells[cix].triangles[t_id].myid = t_id;
                    shells[cix].triangles[t_id].ia = std::stoi(pairs[ 3 ].c_str(), NULL);
                    shells[cix].triangles[t_id].ib = std::stoi(pairs[ 4 ].c_str(), NULL);
                    shells[cix].triangles[t_id].ic = std::stoi(pairs[ 5 ].c_str(), NULL);

                    shells[cix].triangles[t_id].an[0] = strtod(pairs[ 6 ].c_str(), NULL);
                    shells[cix].triangles[t_id].an[1] = strtod(pairs[ 7 ].c_str(), NULL);
                    shells[cix].triangles[t_id].an[2] = strtod(pairs[ 8 ].c_str(), NULL);

                    shells[cix].triangles[t_id].L2[0] = strtod(pairs[ 9 ].c_str(), NULL);
                    shells[cix].triangles[t_id].L2[1] = strtod(pairs[ 10 ].c_str(), NULL);
                    shells[cix].triangles[t_id].L2[2] = strtod(pairs[ 11 ].c_str(), NULL);

                    shells[cix].triangles[t_id].ki[0] = strtod(pairs[ 12 ].c_str(), NULL);
                    shells[cix].triangles[t_id].ki[1] = strtod(pairs[ 13 ].c_str(), NULL);
                    shells[cix].triangles[t_id].ki[2] = strtod(pairs[ 14 ].c_str(), NULL);

                    shells[cix].triangles[t_id].ci[0] = strtod(pairs[ 15 ].c_str(), NULL);
                    shells[cix].triangles[t_id].ci[1] = strtod(pairs[ 16 ].c_str(), NULL);
                    shells[cix].triangles[t_id].ci[2] = strtod(pairs[ 17 ].c_str(), NULL);
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


void Restarter::addBHinges(std::vector<Shell>& shells, int cix) const
{
    std::ifstream os;
    os.open(topologyFile, std::ifstream::in);
    std::string line;

    if ( os.is_open() )
    {
        while ( std::getline (os, line) )
        {
            if ( line.find("SHELLHINGE ") == 0 )
            {
                std::vector<std::string> pairs = split(line, ' ');

                int shell_id = std::stoi(pairs[ 1 ].c_str(), NULL);

                if (shell_id == cix)
                {
                    int b_id = std::stoi(pairs[ 2 ].c_str(), NULL);
                    shells[cix].bhinges[b_id].setId(b_id);

                    shells[cix].bhinges[b_id].D = strtod(pairs[ 3 ].c_str(), NULL);
                    shells[cix].bhinges[b_id].sinTheta0 = strtod(pairs[ 4 ].c_str(), NULL);
                    shells[cix].bhinges[b_id].theta0 = strtod(pairs[ 5 ].c_str(), NULL);
                    shells[cix].bhinges[b_id].x1 = std::stoi(pairs[ 6 ].c_str(), NULL);
                    shells[cix].bhinges[b_id].x2 = std::stoi(pairs[ 7 ].c_str(), NULL);
                    shells[cix].bhinges[b_id].x3 = std::stoi(pairs[ 8 ].c_str(), NULL);
                    shells[cix].bhinges[b_id].x4 = std::stoi(pairs[ 9 ].c_str(), NULL);
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

void Restarter::registerVMap()
{
    std::ifstream os;
    os.open(topologyFile, std::ifstream::in);
    std::string line;

    if ( os.is_open() )
    {
        while ( std::getline (os, line) )
        {
            if ( line.find("VMAP ") == 0 )
            {
                std::vector<std::string> pairs = split(line, ' ');

                int ci = std::stoi(pairs[ 2 ].c_str(), NULL);
                int vi = std::stoi(pairs[ 3 ].c_str(), NULL);
                std::pair<std::string, std::pair<int, int> > new_element(pairs[1], {ci, vi});
                vmap.insert( new_element );

            }
        }
    }
    else
    {
        // print error
    }

    os.close();

}

void Restarter::readLastFrame(std::vector<Shell>& shells) const
{
    std::ifstream os;
    os.open(lastFrameFile, std::ifstream::in);
    std::string line;

    if ( os.is_open() )
    {
        while ( std::getline (os, line) )
        {

            std::vector<std::string> pairs = split(line, ' ');

            if (pairs.size() > 1)
            {
                std::string vkey = pairs[ 0 ].c_str();

                double x = strtod(pairs[ 1 ].c_str(), NULL);
                double y = strtod(pairs[ 2 ].c_str(), NULL);
                double z = strtod(pairs[ 3 ].c_str(), NULL);

                std::pair<int, int> value = vmap.at( vkey );

                int ci = value.first;
                int vi = value.second;

                shells[ci].vertices[vi].r_c.x = x;
                shells[ci].vertices[vi].r_c.y = y;
                shells[ci].vertices[vi].r_c.z = z;
            }
        }

    }
    else
    {
        // print error
    }

    os.close();
}

void Restarter::readFrame(std::string trajFile, std::vector<Shell>& shells, int frameNumber) const
{
    std::ifstream os;
    os.open(trajFile, std::ifstream::in);
    std::string line;

    int frames_counter = 0;

    if ( os.is_open() )
    {
        while ( std::getline (os, line) )
        {

            std::vector<std::string> pairs = split(line, ' ');

            if (pairs.size() == 1)
            {
                frames_counter++;
                continue;
            }

            if (pairs.size() > 1 && frameNumber == frames_counter )
            {
                std::string vkey = pairs[ 0 ].c_str();

                double x = strtod(pairs[ 1 ].c_str(), NULL);
                double y = strtod(pairs[ 2 ].c_str(), NULL);
                double z = strtod(pairs[ 3 ].c_str(), NULL);

                std::pair<int, int> value = vmap.at( vkey );

                int ci = value.first;
                int vi = value.second;

                shells[ci].vertices[vi].r_c.x = x;
                shells[ci].vertices[vi].r_c.y = y;
                shells[ci].vertices[vi].r_c.z = z;
            }
        }
    }
    else
    {
        // print error
    }

    os.close();
}


void Restarter::assignTurgors(std::string turgor_line, std::vector<Shell>& shells) const
{
    std::vector<std::string> pairs = split(turgor_line, ' ');
    double turgor;

    for (uint i = 0; i < shells.size(); i++)
    {
        turgor = strtod(pairs[ 4 + 4 * i + 1 ].c_str(), NULL);
        shells[i].pushDp(turgor);
    }
}

void Restarter::assignBoxSize(std::string box_line, Box& box) const
{
    std::vector<std::string> pairs = split(box_line, ' ');
    double x, y, z;
    x = strtod(pairs[ 0 ].c_str(), NULL);
    y = strtod(pairs[ 1 ].c_str(), NULL);
    z = strtod(pairs[ 2 ].c_str(), NULL);

    box.setX(x);
    box.setY(y);
    box.setZ(z);
}