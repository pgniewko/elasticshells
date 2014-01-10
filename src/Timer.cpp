#include "Timer.h"

using namespace std;

Timer::Timer() 
{
    t0 = 0;
    final = 0;
}

Timer::~Timer()
{
}
    
void Timer::tic() 
{
    t0 = clock();
}

void Timer::toc() 
{
    final += clock() - t0;
}

double Timer::time() 
{
    return double(final) / double(CLOCKS_PER_SEC);
}

inline void print_time() 
{
    time_t rawtime;
    struct tm * timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    cout<<endl<<asctime(timeinfo)<<endl;
}