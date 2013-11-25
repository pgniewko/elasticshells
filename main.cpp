/* 
 * File:   main.cpp
 * Author: Pawel Gniewek
 *
 * Created on November 25, 2013, 11:24 AM
 */
#include <iostream>
#include <stdlib.h>
#include <cstdlib>
#include "./src/Cell.h"


//const char *argp_program_version ="main 0.0.1";

//const char *argp_program_bug_address ="<gniewko.pablo@gmail.com>";

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {
    int N = 5;
    Cell cells[N];
    for (int i = 0; i < N; i++) {
        cells[i].set_coor(1,1,1);
        cells[i].set_radius(1.0);
        cout << cells[i].calc_volume() << endl;
    }
    return 0;
}

