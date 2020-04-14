# Author : Pawel Gniewek (UC Berkeley)
# Email  : pawel.gniewek@berkeley.edu
# License: BSD

# C++ compiler
CXX      := g++-9

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
INCLUDE  := -I/usr/local/include -I$(CURDIR)/src -I/Users/pawel/include
#LIB      := -L/usr/local/lib -L/usr/lib
LIB      := -L/usr/local/lib
DFLAGS   :=

CXXFLAGS := -lm -Wall -O3 -Dulong="unsigned long" -Duint="unsigned int" -std=gnu++11 -fpermissive $(INCLUDE)
LDFLAGS  := $(LIB)
LDLIBS   := -ldl -largp -lsteinhardt -lrndmesh -lnblists -lgfortran -lgsl -lgslcblas


# Local dirs 
SRC      := ./src
BIN      := ./bin

# Installation directory
PREFIX   := /usr/local
