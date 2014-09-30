#ifndef TIMER_H
#define TIMER_H

#include <time.h>
#include <iostream>

#include "Environment.h"
using namespace std;

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
        clock_t t0;
        clock_t final;
};

inline void print_time()
{
    time_t rawtime;
    struct tm* timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    cerr << endl << asctime(timeinfo) << endl;
}

#endif