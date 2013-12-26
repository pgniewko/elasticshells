include config.mk

PROGRAM=$(BIN)/biofilm
TEST_RUNNER=$(TEST_SRC)/test

SOURCES=main.cpp \
        $(wildcard $(SRC)/*.cpp)

OBJECTS=$(SOURCES:.cpp=.o)

TEST_SOURCES=$(wildcard $(TEST_SRC)/*.cpp) \
             $(wildcard $(SRC)/*.cpp)

TEST_OBJECTS=$(TEST_SOURCES:.cpp=.o)

$(PROGRAM): $(OBJECTS)
	g++ -lm $^ -o $@

$(TEST_RUNNER): $(TEST_OBJECTS) 
	g++ -lm $^ $(LDFLAGS) $(LDLIBS) -o $@


main.o: main.cpp

$(SRC)/%.o: $(SRC)/%.cpp $(SRC)/%.h

$(TEST_SRC)/%.o: $(TEST_SRC)/%.cpp $(TEST_SRC)/%.h

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@



# Tell make that these are phony targets
.PHONY: build clean test

#test: clean build
test: $(TEST_RUNNER)
	@$(TEST_RUNNER)
	@echo Test done.

build: $(PROGRAM)
	@echo Build done.

clean:
	rm -f $(PROGRAM) $(TEST_RUNNER) $(OBJECTS) $(TEST_OBJECTS)
	@echo Clean done.