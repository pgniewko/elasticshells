#ifndef TIMER_H
#define TIMER_H

#include <time.h>
#include <sys/time.h>
#include <iostream>

#include "Environment.h"

class Timer
{
    public:
        Timer();
        Timer(const Timer& orig);
        virtual ~Timer();

        void tic();
        void toc();
        double time();

    private:
        clock_t t0 = 0;
        clock_t final = 0;
};

inline void print_time()
{
    time_t rawtime;
    struct tm* timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    std::cerr << std::endl << asctime(timeinfo) << std::endl;
}

double read_timer();

#endif