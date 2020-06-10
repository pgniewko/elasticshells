![Shells](./assets/Studium.jpeg)

>Notice: This is research code that will not necessarily be maintained in the future.
>The code is under development so make sure you are using the most recent version.
>We welcome bug reports and PRs but make no guarantees about fixes or responses.

DESCRIPTION
==================================================
The ```elasticshells``` simulator is written in C/C++ and designed to simulate compacted packings of elastic shells.
For the model description, physics foundations, and implementation details, check the **REFERENCES** section. 

GETTING THE CODE
==================================================
* To get the code:
```
git clone git@github.com:pgniewko/elasticshells.git
```

* To obtain the most recent version of the code:
```
git pull origin master
```

COMPILING AND INSTALLATION - LINUX
==================================================

* Executables
To build the exe file, run in the command line (in the root folder):
```
make
```
To change compilation settings, edit the `config.mk` file.  
The exe file is located in the `./bin` directory.  
The software has been tested with GNU compilers >= 4.8.5.    
To built on OS X systems, the `arpg.h` library must be installed and linked (`LDLIBS   := ... largp ...`) in the `config.mk` file.

* Installation

To install the binary files run in the command line, in the root folder:
```
make install
```

The binary file is installed in the `$PREFIX=/usr/local` (by default) directory.
To uninstall the binary file, run in the command line:
```
make uninstall
```
Please note that you must be root in order to install at `/usr/local`.

* Performance profiling with `gprof`:
    + Compile the code(main.cpp) file with the option: 
        `-pg`
    + Run the binary file as usual: 
        `biofilm [OPTIONS ...]`
    + Run gprof tool in order to analyze the performance:
        `gprof elasticshells gmon.out | less`

* Performance profiling with `valgrind`:
    + Run the binary file with `valgrind --tool=callgrind`:
        `valgrind --tool=callgrind elasticshells`
    + Check the results with the `kcachegrind`:
        ` kcachegrind callgrind.out.xxx`

* Memory check with `valgrind`:
    + Run the binary file:
        `valgrind --leak-check=yes --show-leak-kinds=all biofilm [OPTIONS ...] 2> log`

EXTERNAL LIBRARIES
================
Make sure that the dependencies are compiled with the same compiler (GNU compiler (>=4.8.5) recommended) as `elasticshells` 

* ```argp.h``` library is needed for command line parsing. 
  * On Mac OS X it can be easily installed with:
  ```
  brew install argp-standalone
  ```

  * On linux argp.h is a part of the ```libc6-dev``` library. It can be installed with a command (on Ubuntu):
  ```
  sudo apt-get install libc6-dev
  ```

* GNU Scientific Library ([GSL](http://www.gnu.org/software/gsl/)).
  * On Ubuntu you can install it by running:
  ```
  sudo apt-get install libgsl0-dev
  ```
  * On Mac OS you can install it by running:
  ```
  brew install gsl
  ```

* Linked-lists: [nblists](https://github.com/pgniewko/nblists)
* Mesh generator: [rndmesh](https://github.com/pgniewko/rndmesh)
* Order metrics: [steinhardt](https://github.com/pgniewko/steinhardt)
* Fast-math library: [fastmah](https://github.com/pgniewko/fastmath)

USAGE
=====
The code can be used to generate and analyze packings of elastic shells with or without periodic boundary conditions.
Please refer to the ```examples``` directory for scripts to execute the code. 

* ```run_fem.sh``` - Script to compress 5 elastic shells in a rigid box. Change the ```$PBC``` variable for the periodic boundary condition.
* ```restart.sh``` - Resart the simulation setup from the ```run_fem.sh``` script.
* ```analyze.sh``` - Read the trajectory generated by ```run_fem.sh``` script, calculate observables defined in the ```observers.config``` configuration file, and save in an output file.


LICENSE
=======
The library is open-source. If you want to cite the library in any published work please contact me at gniewko.pablo@gmail.com for an information about credits.

COPYRIGHT NOTICE
================
Copyright (C) 2014-2018, Pawel Gniewek  
Email : gniewko.pablo@gmail.com  
All rights reserved.  
License: BSD 3  

REFERENCES
==========
1. "Mechanics of Confined Microbial Populations", P. Gniewek, Ph.D. Thesis, UC Berkeley (2018)


TODO
=====
1. Surface and volume evaluation with Loop shell subdivision method [paper 1](https://github.com/pgniewko/elasticshells/blob/master/assets/peps/COS.2000.pdf) [paper 2](https://github.com/pgniewko/elasticshells/blob/master/assets/peps/FK.2006.pdf)       
2. Add thermal fluctuations and negative internal pressures (to stuck shells buclking) [paper](https://github.com/pgniewko/elasticshells/blob/master/assets/peps/SKB.2020.pdf) [SI](https://github.com/pgniewko/elasticshells/blob/master/assets/peps/SKB.2020-SI.pdf)          
3. Consider also Oriented Particle Systems ([OPS](https://github.com/pgniewko/elasticshells/blob/master/assets/peps/ST.1992-OPS.pdf)) model         
4. We could try to simulate this [experiment](https://github.com/pgniewko/elasticshells/blob/master/assets/peps/TDS.2020.pdf).         

ACKNOWLEDGMENTS
===============
We thank Dominik Gront (University of Warsaw) for sharing [BioShell](http://bioshell.pl/) `Logger` and
`LogManager` classes.
