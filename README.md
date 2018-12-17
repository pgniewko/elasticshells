>Notice: This is research code that will not necessarily be maintained in the future.
>The code is under development so make sure you are using the most recent version.
>We welcome bug reports and PRs but make no guarantee about fixes or responses.

DESCRIPTION
==================================================
```elasticshells``` simulator is written in C/C++, designed to simulate compacted packings of elastic shells.

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
To build the exe file run in the command line:
```
make
```
To change compilation settings edit `config.mk` file.  
Exe file is located in `./bin` directory.  
The software has been tested with GNU compilers >= 4.8.5.    
To built on OS X systems `arpg.h` library must be installed and linked (`LDLIBS   := ... largp ...`) in `config.mk` file.

* Installation

To install binary files run in the command line:
```
make install
```

The binary file is installed in `$PREFIX=/usr/local` (by default) directory.
To uninstall the binary file run in the command line:
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
Make sure that the dependencies are compiled with the same compiled (GNU compiler (>=4.8.5) recommended) as `elasticshells` 

* ```argp.h``` library is needed for command line parsing. 
  * On Mac OS X it can be easily installed with:
  ```
  brew install argp-standalone
  ```

  * On linux argp.h is a part of ```libc6-dev``` library. It can be installed with a command (on Ubuntu):
  ```
  sudo apt-get install libc6-dev
  ```

* GNU Scientific Library ([GSL](http://www.gnu.org/software/gsl/)).
  * On Ubuntu you can install it by running:
  ```
  apt-get install libgsl0-dev
  ```
  * On Mac OS you can install it by runnig:
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

* ```run_fem.sh``` - 
* ```restart.sh``` - 
* ```analyze.sh``` - 


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

ACKNOWLEDGMENTS
===============
We thank Dominik Gront (University of Warsaw) for sharing [BioShell](http://bioshell.pl/) `Logger` and
`LogManager` classes.
