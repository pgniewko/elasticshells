include config.mk

PROGRAM      := $(BIN)/biofilm
TEST_RUNNER  := $(TESTS)/testsrunner

SOURCES      := main.cpp \
                $(wildcard $(SRC)/*.cpp) \
                $(wildcard $(SRC)/geometry/*.cpp) \
                $(wildcard $(SRC)/simulation/*.cpp) \
		$(wildcard $(SRC)/geometry/algorithms/*.cpp) \
		$(wildcard $(SRC)/force/*.cpp) \
		$(wildcard $(SRC)/utils/io/*.cpp) \

HEADERS      := $(wildcard $(SRC)/*.h) \
		$(wildcard $(SRC)/exceptions/*.h) \
	        $(wildcard $(SRC)/geometry/*.h) \
	        $(wildcard $(SRC)/simulation/*.h) \
	        $(wildcard $(SRC)/geometry/algorithms/*.h) \
		$(wildcard $(SRC)/force/*.h) \
	        $(wildcard $(SRC)/utils/io/*.h) \

TEST_SOURCES := $(wildcard $(TESTS)/*.cpp) \
                $(wildcard $(SRC)/*.cpp) \
		$(wildcard $(SRC)/exceptions/*.cpp) \
                $(wildcard $(SRC)/geometry/*.cpp) \
                $(wildcard $(SRC)/simulation/*.cpp) \
		$(wildcard $(SRC)/geometry/algorithms/*.cpp) \
		$(wildcard $(SRC)/force/*.cpp) \


OBJECTS      := $(SOURCES:.cpp=.o)


TEST_OBJECTS := $(TEST_SOURCES:.cpp=.o)

#Linking commands:
$(PROGRAM): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@

$(TEST_RUNNER): $(TEST_OBJECTS) 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)

#Compilation commands:
main.o: $(SOURCES) $(HEADERS) #main.cpp 

$(SRC)/%.o: $(SRC)/%.cpp $(SRC)/%.h

$(TESTS)/%.o: $(TESTS)/%.cpp $(TESTS)/%.h

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@


# Tell make that these are phony targets
.PHONY: build clean test


build: $(PROGRAM)
	@echo Build done.

install:
	@echo You must be root to install.
	
clean:
	rm -f $(PROGRAM) $(TEST_RUNNER) $(OBJECTS) $(TEST_OBJECTS)
	@echo Clean done.
	
tests: $(TEST_RUNNER)
	@$(TEST_RUNNER)
	@echo Test done.

indent:
	@astyle --style=allman -r -xl -C -xG -SKNL -wfpHj -k1 "*.cpp"
	@astyle --style=allman -r -xl -C -xG -SKNL -wfpHj -k1 "*.h"
	@./astyle-clean.sh
	
