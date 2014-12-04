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

RELATED LIBRARIES/PROGRAMS
================
* For compiling and running tests CppUnit [link](sourceforge.net/projects/cppunit) is needed.

* GNU Scientific Library (GSL) [link](http://www.gnu.org/software/gsl/)
On Ubuntu you can install it by running:
```
apt-get install libgsl0-dev
```

* Steinhardt order parameters library [link](https://github.com/nquesada/steinhardt).



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
