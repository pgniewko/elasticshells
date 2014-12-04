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

RELATED PROGRAMS
================
Biofilm software is using external library(steinhardt [link](https://github.com/nquesada/steinhardt)) to calculate order parameters.
In order to make it running you need GNU Scientific Library (GSL).
On Ubuntu you can install it with the commad:
```
apt-get install libgsl0-dev
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
