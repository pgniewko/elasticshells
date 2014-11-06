Description
==================================================

Compilation - Linux / Mac OS / Windows 
==================================================

* Performance profiling with `gprof`:
    + Compile the code(main.cpp) file witg the option: 
        `-pg`
    + Run the binary file as usual: 
        `biofilm [OPTIONS ...]`
    + Run gprof tool in order to analyze the performance:
        `gprof ./bin/biofilm gmon.out | less`

Related programs
================


Contents
========

* To check the number of lines of code run
```
    find . -name "*.cpp" -exec wc {} \; | awk 'BEGIN{SUM=0} {SUM += $1} END{print SUM}'
    find . -name "*.h"   -exec wc {} \; | awk 'BEGIN{SUM=0} {SUM += $1} END{print SUM}'
```


Usage
=====


Copyright Notice
================
Copyright (C) 2014 -,  Pawel Gniewek
All rights reserved.
License: BSD

Acknowledgments
===============
