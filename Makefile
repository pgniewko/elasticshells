include config.mk

PROGRAM      := $(BIN)/biofilm
TEST_RUNNER  := $(TESTS)/test

SOURCES      := main.cpp \
                $(wildcard $(SRC)/*.cpp) \
                $(wildcard $(SRC)/geometry/*.cpp) \
		$(wildcard $(SRC)/geometry/algorithms/*.cpp) \
		$(wildcard $(SRC)/force/*.cpp) \

HEADERS      := $(wildcard $(SRC)/*.h) \
	        $(wildcard $(SRC)/geometry/*.h) \
	        $(wildcard $(SRC)/geometry/algorithms/*.h) \
		$(wildcard $(SRC)/force/*.h) \

OBJECTS      := $(SOURCES:.cpp=.o)

TEST_SOURCES := $(wildcard $(TESTS)/*.cpp) \
                $(wildcard $(SRC)/*.cpp)

TEST_OBJECTS := $(TEST_SOURCES:.cpp=.o)

#Linking commands:
$(PROGRAM): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@

$(TEST_RUNNER): $(TEST_OBJECTS) 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LDLIBS) $^ -o $@

#Compilation commands:
main.o: $(SOURCES) $(HEADERS) #main.cpp 

$(SRC)/%.o: $(SRC)/%.cpp $(SRC)/%.h

$(TESTS)/%.o: $(TESTS)/%.cpp $(TESTS)/%.h

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


# Tell make that these are phony targets
.PHONY: build clean test

test: $(TEST_RUNNER)
	@$(TEST_RUNNER)
	@echo Test done.

build: $(PROGRAM)
	@echo Build done.

clean:
	rm -f $(PROGRAM) $(TEST_RUNNER) $(OBJECTS) $(TEST_OBJECTS)
	@echo Clean done.
	
indent:
	@astyle --style=allman -r -xl -C -xG -SKNL -wfpHj -k1 "*.cpp"
	@astyle --style=allman -r -xl -C -xG -SKNL -wfpHj -k1 "*.h"
	@./astyle-clean.sh
	
install:
	@echo You must be root to install.