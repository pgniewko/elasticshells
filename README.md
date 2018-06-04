DESCRIPTION
==================================================

Biofilm simulator is written in C/C++ and uses modern computer
architectures and technologies: OpenMP for shared-memory systems, 
SSE vectorization for x86_64 CPU.

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
To build executables run:
```
make
```
To change compilation settings edit `config.mk` file.  
Executables will be located in `./bin` directory.  
For a correct compilation g++ >= 4.9 is necessary.  
To built on OS X systems `arpg.h` library must be installed and linked (`LDLIBS   := ... largp ...`) in `config.mk` file.

* Installation

To install binary files run:
```
make install
```

Binary files will be installed in `$PREFIX=/usr/local` (by default) directory.
To uninstall executables use:
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
        `gprof biofilm gmon.out | less`

* Performance profiling with `valgrind`:
    + Run the binary file with `valgrind --tool=callgrind`:
        `valgrind --tool=callgrind biofilm`
    + Check the results with the `kcachegrind`:
        ` kcachegrind callgrind.out.xxx`

* Memory check with `valgrind`:
    + Run the binary file:
        `valgrind --leak-check=yes --show-leak-kinds=all biofilm [OPTIONS ...] 2> log`

EXTERNAL LIBRARIES
================
* For compiling and running tests [CppUnit](sourceforge.net/projects/cppunit) is needed.
* ```argp.h``` library is needed for command line parsing. On Mac OS X it can be easily installed with:
```
brew install argp-standalone
```

On linux argp.h is a part of ```libc6-dev``` library. It can be installed with a command (on Ubuntu):
```
sudo apt-get install libc6-dev
```

* GNU Scientific Library ([GSL](http://www.gnu.org/software/gsl/)).
On Ubuntu you can install it by running:
```
apt-get install libgsl0-dev
```

* Steinhardt library for order parameters parameters calculations.
The original code can be downloaded from [github](https://github.com/nquesada/steinhardt),
or from the local [copy](https://bitbucket.org/pawelgniewek/steinhardt). 

* Linked-lists: [nblists](https://github.com/pgniewko/nblists)
* Mesh generator: [rndmesh](https://github.com/pgniewko/rndmesh)
* Order metrics: [steinhardt](https://github.com/pgniewko/steinhardt)
* Fast-math library: [fastmah](https://github.com/pgniewko/fastmath)

CONTENTS
========

* To check the number of lines of code run
```
find . -name "*.cpp" -exec wc {} \; | awk 'BEGIN{SUM=0} {SUM += $1} END{print SUM}'
find . -name "*.h"   -exec wc {} \; | awk 'BEGIN{SUM=0} {SUM += $1} END{print SUM}'
```


USAGE
=====


LICENSE
=======
The library is open-source. If you want to cite the library in any published work please contact me at gniewko.pablo@gmail.com 
for an information about credits or check my [website](http://meetpg.pl). Physics background and models benchmarking can be fount in the [notes](http://meetpg.pl/notes.html).

COPYRIGHT NOTICE
================
Copyright (C) 2014-2018,  Pawel Gniewek  
Email  : gniewko.pablo@gmail.com  
All rights reserved.  
License: BSD 3  

ACKNOWLEDGMENTS
===============
We thank Dominik Gront (University of Warsaw) for sharing [BioShell](http://bioshell.pl/) `Logger` and
`LogManager` classes.
