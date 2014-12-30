include config.mk

TARGET       := $(BIN)/biofilm
TEST_RUNNER  := $(TESTS)/testsrunner


SOURCES	     := main.cpp \
		$(shell find $(SRC) -type f -name "*.cpp")


HEADERS	     := $(shell find $(SRC) -type f -name "*.h")


TEST_SOURCES := $(shell find $(TESTS) -type f -name "*.cpp") \
		$(shell find $(SRC)   -type f -name "*.cpp")


OBJECTS      := $(SOURCES:.cpp=.o)


TEST_OBJECTS := $(TEST_SOURCES:.cpp=.o)

#DEPS         := $(SOURCES:.cpp=.d)
DEPS         := $(OBJECTS:.o=.d)

#Linking commands:
$(TARGET): $(OBJECTS)
	@echo Linking...
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)

$(TEST_RUNNER): $(TEST_OBJECTS) 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)

#Compilation commands:
main.o: $(SOURCES) $(HEADERS) #main.cpp 

$(SRC)/%.o: $(SRC)/%.cpp $(SRC)/%.h

$(TESTS)/%.o: $(TESTS)/%.cpp $(TESTS)/%.h

%.o: %.cpp %.h
	$(CXX) -MMD -MP $(CXXFLAGS) -c $< -o $@
	
# Tell make that these are phony targets
.PHONY: build clean tests install uninstall indent


build: $(TARGET)
	@echo Build done.

clean:
	@echo Cleaning...
	rm -f $(TARGET) $(TEST_RUNNER) $(OBJECTS) $(TEST_OBJECTS) $(DEPS)
	@echo Clean done.
	
tests: $(TEST_RUNNER)
	@$(TEST_RUNNER)
	@echo Test done.

install: $(TARGET)
	@echo You must be root to install. Have password ready!
	sudo install -m 755 $(TARGET) $(PREFIX)/bin
	@echo "Installation complete!"

uninstall:
	@echo You must be root to uninstall. Have password ready!
	@if [ -f $(PREFIX)/bin/biofilm ]; \
	then \
	    sudo rm $(PREFIX)/bin/biofilm; \
	fi
	@echo '$(PREFIX)/bin/biofilm' uninstalled!
	
indent:
	@astyle --style=allman -r -xl -C -xG -SKNL -wfpHj -k1 "*.cpp"
	@astyle --style=allman -r -xl -C -xG -SKNL -wfpHj -k1 "*.h"
	@./astyle-clean.sh

-include $(DEPS)