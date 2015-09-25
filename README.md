DESCRIPTION
Biofilm simulator is written in C/C++ and uses modern computer
architectures and technologies: OpenMP for shared-memory systems, 
SSE vectorization for x86_64 CPU.
==================================================
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
To built on OS X systems `arpg.h` library must be installed and
linked (`LDLIBS   := ... largp ...`) in `config.mk` file.

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
        `gprof ./bin/biofilm gmon.out | less`

* Performance profiling with `valgrind`:
    + Run the binary file with `valgrind --tool=callgrind`:
        `valgrind --tool=callgrind ./bin/biofilm`
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


* GNU Scientific Library ([GSL](http://www.gnu.org/software/gsl/)).
On Ubuntu you can install it by running:
```
apt-get install libgsl0-dev
```

* Steinhardt order parameters [library](https://github.com/nquesada/steinhardt).
The library can be found in `lib/steinhardt` directory.
First we compile a static library, then combine objects into library, create an index within the library,
and finally we can install it in `/usr/local/lib`:
```
$ g++ -O3 -Wall -lm -lgsl -lgslcblas -c -o steinhardt.o steinhardt.c 
$ ar -rcs libsteinhardt.a steinhardt.o
$ ranlib libsteinhardt.a
$ sudo cp libsteinhardt.a /usr/local/lib
$ sudo cp steinhardt.h /usr/local/include
```


CONTENTS
========

* To check the number of lines of code run
```
find . -name "*.cpp" -exec wc {} \; | awk 'BEGIN{SUM=0} {SUM += $1} END{print SUM}'
find . -name "*.h"   -exec wc {} \; | awk 'BEGIN{SUM=0} {SUM += $1} END{print SUM}'
```


USAGE
=====


COPYRIGHT NOTICE
================
Copyright (C) 2014-2015,  Pawel Gniewek
All rights reserved.
License: BSD

ACKNOWLEDGMENTS
===============
We thank Dominik Gront (University of Warsaw) for sharing BioShell `Logger` and
`LogManager` classes.
