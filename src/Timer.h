#ifndef TIMER_H
#define TIMER_H

#include <time.h>
#include <iostream>

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
	
#endif

