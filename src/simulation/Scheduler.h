#ifndef SCHEDULER_H
#define	SCHEDULER_H

#include <iostream>
#include <vector>
#include <sstream>
#include <string>           // std::stoi
#include <climits>          // MAX_INT
#include <cmath>            // std::max

#include "utils/Logger.h"
#include "utils/utils.h"

struct schedule_t
{
    int n_steps;
    int interval;
    double dx;
    double dy;
    double dz;
    double rx;
    double ry;
    double rz;
    int counter;
    double vf;

    schedule_t() : n_steps(0), interval(1), dx(0), dy(0), dz(0), rx(0), ry(0), rz(0), counter(0), vf(0) {}

    schedule_t(int ns, int in, double _dx, double _dy, double _dz, double _rx, double _ry, double _rz, double _vf) :
        n_steps(ns), interval(in), dx(_dx), dy(_dy), dz(_dz), rx(_rx), ry(_ry), rz(_rz), counter(0), vf(_vf) {}

};

class Scheduler
{
    public:
        Scheduler();
        Scheduler(const Scheduler& orig);
        virtual ~Scheduler();
        std::vector<std::string> read_schedule_file();
        void set_file_name(std::string);
        void configure_schedule();
        void register_schedules();
        void check_schedule();
        void change_schedule();
        void print_schedule();
        void save_remaining_schedule();
        void set_default(int, int, double, double, double, double, double, double);
        void execute(double&, double&, double&, const double);

        bool nth_todo();

    protected:
        std::string schedulefile;
        std::vector<schedule_t> schedules;
        schedule_t default_schedule;
        schedule_t current_schedule;

    private:
        bool schedule_registered = false;
        static utils::Logger schedule_logger;

};

#endif	/* SCHEDULER_H */

