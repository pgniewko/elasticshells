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
INCLUDE  := -Isrc
LIB      := -Llib

CXXFLAGS := -lm -Wall -O3 -fpermissive $(DEBUG) $(INCLUDE) -Wunused-but-set-variable
LDFLAGS  := $(LIB)
LDLIBS   := -lcppunit -ldl

# Local dirs 
TESTS    := ./tests
SRC      := ./src
BIN      := ./bin

# Installation directory
#PREFIX=/usr/local
