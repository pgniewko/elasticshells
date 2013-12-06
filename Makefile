# Author : Pawel Gniewek (UC Berkeley)
# Email  : pawel.gniewek@berkeley.edu
# License: BSD

include config.mk

PROGRAM=$(BIN)/biofilm
SOURCES=./main.cpp \
  $(SRC)/Cell.cpp \
  $(SRC)/YeastCell.cpp \
  $(SRC)/Simulator.cpp

OBJECTS=$(SOURCES:.cpp=.o)

$(PROGRAM): $(OBJECTS)
	g++ -lm $^ -o $@

main.o: main.cpp

Cell.o: $(SRC)/Cell.cpp $(SRC)/Cell.h

Simulator.o: $(SRC)/Simulator.cpp $(SRC)/Simulator.h

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ 


# Tell make that these are phony targets
.PHONY: build clean test

test: clean build
	@echo Test done.

build: $(PROGRAM)
	@echo Build done.

clean:
	rm -f $(PROGRAM) $(OBJECTS)
	@echo Clean done.