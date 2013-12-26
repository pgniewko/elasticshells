# Author : Pawel Gniewek (UC Berkeley)
# Email  : pawel.gniewek@berkeley.edu
# License: BSD

# This a common configuration file that includes definitions used by all
# the Makefiles.

# C++ compiler
CXX=g++

# Flags for the C++ compiler
# The option "-Wall" tells the compiler to print all warnings.
# "-O3" - optimization level
DEBUG=-g
CXXFLAGS=-Wall -O3 -fpermissive $(DEBUG)
LDLIBS=-lcppunit

# Local dirs 
SRC=./src
LIB=./lib
BIN=./bin

TEST_SRC=./src/test

# Relative include and library paths for compilation of the examples
#E_INC=-I.
#E_LIB=-L.

# Installation directory
#PREFIX=/usr/local

# Install command
#INSTALL=install

# Flags for install command for executable
#IFLAGS_EXEC=-m 0755

# Flags for install command for non-executable files
#IFLAGS=-m 0644
