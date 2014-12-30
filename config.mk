# Author : Pawel Gniewek (UC Berkeley)
# Email  : pawel.gniewek@berkeley.edu
# License: BSD

# C++ compiler
CXX      := g++

# Flags for the C++ compiler
# -g --Enable debugging
# -Wall --Turn on all warnings
# -D_USE_FIXED_PROTOTYPES_
# --Force the compiler to use the correct headers
# -ansi --Don't use GNU ext; do use ansi standard.
# "-O3" - optimization level
DEBUG    := -g

# Relative include and library paths for compilation of the examples
INCLUDE  := -I/usr/local/include -I$(CURDIR)/src
LIB      := -L/usr/lib -L/usr/local/lib
DFLAGS   := #-DTESTS

CXXFLAGS := -lm -Wall -O3 -std=gnu++0x -fpermissive $(DEBUG) $(INCLUDE) $(DFLAGS)
LDFLAGS  := $(LIB)
LDLIBS   := -lcppunit -ldl -lsteinhardt -lgsl -lgslcblas


# Local dirs 
TESTS    := ./tests
SRC      := ./src
BIN      := ./bin

# Installation directory
PREFIX   := /usr/local
