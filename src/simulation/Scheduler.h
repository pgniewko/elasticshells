#ifndef SCHEDULER_H
#define	SCHEDULER_H

#include <iostream>
#include <vector>
#include <sstream>
#include <string>           // std::stoi

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
    
    schedule_t() : n_steps(0), interval(0), dx(0), dy(0), dz(0), rx(0), ry(0), rz(0), counter(0) {}
    
    schedule_t(int ns, int in, double _dx, double _dy, double _dz, double _rx, double _ry, double _rz) :
    n_steps(ns), interval(in), dx(_dx), dy(_dy), dz(_dz), rx(_rx), ry(_ry), rz(_rz), counter(0) {}
    
};

class Scheduler
{
    public:
        Scheduler();
        //Scheduler(char*);
        Scheduler(const Scheduler& orig);
        virtual ~Scheduler();
        std::vector<std::string> readScheduleFile();
        void setFileName(char*);
        void configureSchedule();
        void registerSchedules();
        void checkSchedule();
        void changeSchedule();
        void printSchedule();
        void setDefault(int, int, double, double, double, double, double, double);
        void execute(double&, double&, double&);
        //int total_time;
        //int recent_time;
                   
    protected:
        std::string schedulefile;
        std::vector<schedule_t> schedules;
        schedule_t default_schedule;
        schedule_t current_schedule;
        
private:
    static utils::Logger schedule_logger;    

};

#endif	/* SCHEDULER_H */

