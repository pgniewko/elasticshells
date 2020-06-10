#include "Scheduler.h"

utils::Logger Scheduler::schedule_logger("schedule_logger");

Scheduler::Scheduler() {}

Scheduler::Scheduler(const Scheduler& orig) :
    schedulefile(orig.schedulefile), default_schedule(orig.default_schedule), current_schedule(orig.current_schedule)
{
    for (uint i = 0; i < orig.schedules.size(); i++)
    {
        schedules.push_back(orig.schedules[i]);
    }
}

Scheduler::~Scheduler() {}

std::vector<std::string> Scheduler::read_schedule_file()
{
    std::ifstream cfile;
    cfile.open(schedulefile, std::ifstream::in);
    std::vector<std::string> list;
    std::string line;

    if ( cfile.is_open() )
    {
        while ( std::getline (cfile, line) )
        {
            if ( !(line.at(0) == '#') && !(line.at(0) == ' ') )
            {
                list.push_back(line);
            }
        }
    }
    else
    {
        schedule_logger << utils::LogLevel::WARNING << "Observers configuration file COULD NOT BE FOUND" << "\n";
    }

    cfile.close();
    return list;
}

void Scheduler::register_schedules()
{
    std::vector<std::string> list = read_schedule_file();

    if ( list.size() > 0)
    {
        schedule_registered = true;
    }

    std::vector<std::string> single_line;

    for (std::vector<std::string>::iterator it = list.begin(); it != list.end(); ++it)
    {
        single_line = split( *it, ' ');

        if (single_line.size() >= 8)
        {

            int ns = std::stoi(single_line[ 0 ].c_str(), NULL);
            int i  = std::stoi(single_line[ 1 ].c_str(), NULL);
            double v1 = strtod(single_line[ 2 ].c_str(), NULL);
            double v2 = strtod(single_line[ 3 ].c_str(), NULL);
            double v3 = strtod(single_line[ 4 ].c_str(), NULL);
            double v4 = strtod(single_line[ 5 ].c_str(), NULL);
            double v5 = strtod(single_line[ 6 ].c_str(), NULL);
            double v6 = strtod(single_line[ 7 ].c_str(), NULL);
            double v7 = strtod(single_line[ 8 ].c_str(), NULL);
            v7 = fmin(1.0, fabs( v7 ) );
            schedule_t new_schedule(ns, i, v1, v2, v3, v4, v5, v6, v7);
            schedules.push_back(new_schedule);
        }
    }
}

void Scheduler::set_file_name(std::string schf)
{
    schedulefile = schf;
}

void Scheduler::configure_schedule()
{
    if ( !schedules.empty() )
    {
        schedule_logger << utils::LogLevel::FINEST << "FIRST SCHEDULE ON THE STACK !\n";
        current_schedule = schedules[0];
    }
}

void Scheduler::print_schedule()
{
    for (uint i = 0; i < schedules.size(); i++)
    {
        schedule_t s = schedules[i];
        std::cout << s.n_steps << " " << s.interval << " " << s.dx << " " << s.dy << " " << s.dz << " " << s.rx << " " << s.ry << " " << s.rz << " " << s.vf << std::endl;
    }

    std::cout << "DEFAULT" << std::endl;
    std::cout << default_schedule.n_steps << " " << default_schedule.interval << " " << default_schedule.dx << " " << default_schedule.dy << " " << default_schedule.dz << " " << default_schedule.rx << " " << default_schedule.ry << " " << default_schedule.rz << std::endl;
}

void Scheduler::save_remaining_schedule()
{
    std::ofstream ofile;
    ofile.open(std::string(schedulefile) + std::string(".remain"), std::ifstream::out);

    if ( ofile.is_open() )
    {
        for (uint i = 0; i < schedules.size(); i++)
        {
            schedule_t s = schedules[i];

            if (i == 0)
            {
                ofile << (s.n_steps - current_schedule.counter) << " " << s.interval << " ";
            }
            else
            {
                ofile << s.n_steps << " " << s.interval << " ";
            }

            ofile << s.dx << " " << s.dy << " " << s.dz << " ";
            ofile << s.rx << " " << s.ry << " " << s.rz << " " << s.vf << std::endl;
        }

        ofile.close();

    }
    else
    {
        schedule_logger << utils::LogLevel::WARNING << "Observers' configuration file COULD NOT BE FOUND" << "\n";
    }
}

void Scheduler::set_default(int ns, int in, double _dx, double _dy, double _dz, double _rx, double _ry, double _rz)
{
    default_schedule = schedule_t(ns, in, _dx, _dy, _dz, _rx, _ry, _rz, 0.95);
}

void Scheduler::execute(double& dx, double& dy, double& dz, const double vf_)
{
    if ( !schedules.empty() )
    {
        if (current_schedule.counter % current_schedule.interval == 0)
        {
            dx = current_schedule.dx + uniform(-current_schedule.rx, current_schedule.rx);
            dy = current_schedule.dy + uniform(-current_schedule.ry, current_schedule.ry);
            dz = current_schedule.dz + uniform(-current_schedule.rz, current_schedule.rz);
        }


        if (vf_ >= current_schedule.vf)
        {
            dx = 0.0;
            dy = 0.0;
            dz = 0.0;
        }

        current_schedule.counter++;

        if (current_schedule.counter > current_schedule.n_steps || vf_ >= current_schedule.vf)
        {
            schedules.erase(schedules.begin());

            if ( !schedules.empty() )
            {
                current_schedule = schedules[0];
                current_schedule.counter = 0;
                schedule_logger << utils::LogLevel::FINEST << "NEW SCHEDULE ON THE STACK\n";
            }
            else
            {
                schedule_logger << utils::LogLevel::FINEST << "SCHEDULES LIST REACHED THE END. SWITCH TO THE DEFAULT SCHEDULE\n";
            }
        }
    }
    else
    {
        if (default_schedule.counter == 0)
        {
            schedule_logger << utils::LogLevel::INFO << "NO MORE SCHEDULES ON THE STACK. DEFAULT SCHEDULE IS EXECUTED.\n";
        }

        if ( default_schedule.counter % default_schedule.interval == 0)
        {
            dx = default_schedule.dx + uniform(-default_schedule.rx, default_schedule.rx);
            dy = default_schedule.dy + uniform(-default_schedule.ry, default_schedule.ry);
            dz = default_schedule.dz + uniform(-default_schedule.rz, default_schedule.rz);
        }

        if (vf_ >= default_schedule.vf)
        {
            dx = 0.0;
            dy = 0.0;
            dz = 0.0;
            default_schedule.counter = default_schedule.n_steps; // THIS WILL TERMINATE THE SIMULATION
        }

        default_schedule.counter++;
    }

}

bool Scheduler::nth_todo()
{
    if ( schedules.empty() && schedule_registered )
    {
        return true;
    }

    if ( default_schedule.dx == 0.0 && default_schedule.dy == 0.0 && default_schedule.dz == 0.0 && schedules.empty() )
    {
        return true;
    }

    if ( default_schedule.dx == 0.0 && default_schedule.dy == 0.0 && default_schedule.dz == 0.0 && !schedule_registered )
    {
        return true;
    }

    if ( default_schedule.n_steps == 0 && !schedule_registered )
    {
        return true;
    }

    if ( default_schedule.counter > default_schedule.n_steps )
    {
        return true;
    }

    return false;
}