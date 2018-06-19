# Author : Pawel Gniewek (UC Berkeley)
# Email  : pawel.gniewek@berkeley.edu
# License: BSD

# C compiler
CC	     := gcc

# C++ compiler
CXX      := g++-4.9
# Flags for the C++ compiler
# -g --Enable debugging
# -Wall --Turn on all warnings
# -D_USE_FIXED_PROTOTYPES_
# --Force the compiler to use the correct headers
# -ansi: In C++ mode, it is equivalent to -std=c++98
# "-O3" - optimization level
# -march=native - optimize for the particular architecture
DBGFLAGS := -DDEBUG -g

# Relative include and library paths for compilation of the examples
INCLUDE  := -I/usr/local/include -I$(CURDIR)/src #-I/Users/pawel/include
#LIB      := -L/usr/local/lib -L/usr/lib
LIB      := -L/usr/local/lib
DFLAGS   := -DTESTS

CXXFLAGS := -lm -Wall -O3 -std=gnu++11 -DFAST -Dulong="unsigned long" -Duint="unsigned int" -fopenmp -fpermissive $(INCLUDE)
LDFLAGS  := $(LIB)
LDLIBS   := -lcppunit -ldl -largp -lsteinhardt -lrndmesh -lnblists -lgfortran -lgsl -lgslcblas


# Local dirs 
TESTS    := ./tests
SRC      := ./src
BIN      := ./bin

# Installation directory
PREFIX   := /usr/local
