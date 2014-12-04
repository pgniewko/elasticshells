DESCRIPTION
==================================================
COMPILATION - LINUX / MAX OS / WINDOWS 
==================================================

* Performance profiling with `gprof`:
    + Compile the code(main.cpp) file witg the option: 
        `-pg`
    + Run the binary file as usual: 
        `biofilm [OPTIONS ...]`
    + Run gprof tool in order to analyze the performance:
        `gprof ./bin/biofilm gmon.out | less`

EXTERNAL LIBRARIES
================
* For compiling and running tests [CppUnit](sourceforge.net/projects/cppunit) is needed.

* GNU Scientific Library ([GSL](http://www.gnu.org/software/gsl/)).
On Ubuntu you can install it by running:
```
apt-get install libgsl0-dev
```

* Steinhardt order parameters [library](https://github.com/nquesada/steinhardt).
The library can be found in `lib/steinhardt` directory.
First we need to compile a static library, and then we can install it in `/usr/local/lib`:
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
Copyright (C) 2014 -,  Pawel Gniewek
All rights reserved.
License: BSD

ACKNOWLEDGMENTS
===============
We thank Dominik Gront (Univeristy of Warsaw) for sharing BioShell `Logger` and
`LogManager` classes.