include config.mk

RM          := rm -f
MKDIR	    := mkdir -p

DEBUG ?= 0
ifeq ($(DEBUG), 1)
    CXXFLAGS := $(CXXFLAGS) $(DBGFLAGS)	
else
ifneq ($(DEBUG), 0)
    CXXFLAGS := $(CXXFLAGS) $(DEBUG)	
endif
endif

TEST ?= 0
ifeq ($(TEST), 1)
    CXXFLAGS := $(CXXFLAGS) $(DFLAGS)
endif

TARGET       := $(BIN)/elasticshells
TEST_RUNNER  := $(TESTS)/testsrunner


SOURCES	     := main.cpp \
		$(shell find $(SRC) -type f -name "*.cpp")

HEADERS	     := $(shell find $(SRC) -type f -name "*.h")

TEST_SOURCES := $(shell find $(TESTS) -type f -name "*.cpp") \
		$(shell find $(SRC)   -type f -name "*.cpp")

OBJECTS      := $(SOURCES:.cpp=.o)

TEST_OBJECTS := $(TEST_SOURCES:.cpp=.o)

DEPS         := $(OBJECTS:.o=.d)

#Linking commands:
$(TARGET): $(OBJECTS)
	@echo LINKING ...
	$(MKDIR) $(@D)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)
	@echo BUILDING IS DONE

$(TEST_RUNNER): $(TEST_OBJECTS) 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@ -lcppunit $(LDLIBS)

#Compilation commands:
main.o: $(SOURCES) $(HEADERS) #main.cpp 

$(SRC)/%.o: $(SRC)/%.cpp $(SRC)/%.h

$(TESTS)/%.o: $(TESTS)/%.cpp $(TESTS)/%.h

%.o: %.cpp %.h
	$(CXX) -MMD -MP $(CXXFLAGS) -c $< -o $@
	
# Tell make that these are phony targets
.PHONY: build clean tests install uninstall indent


build: $(TARGET)
	@echo BUILDING IS DONE

clean:
	@echo Cleaning...
	$(RM) $(TARGET) $(TEST_RUNNER) $(OBJECTS) $(TEST_OBJECTS) $(DEPS)
	
tests: $(TEST_RUNNER)
	@$(TEST_RUNNER)
	@echo Test done.

install: $(TARGET)
	@echo You must be root to install. Have password ready!
	sudo install -m 755 $(TARGET) $(PREFIX)/bin
	@echo 'elasticshells has been installed at $(PREFIX)/bin'
	@echo "INSTALLATION COMPLETE !"

uninstall:
	@echo You must be root to uninstall. Have password ready!
	@if [ -f $(PREFIX)/bin/elasticshells ]; \
	then \
	    sudo rm $(PREFIX)/bin/elasticshells; \
	fi
	@echo '$(PREFIX)/bin/elasticshells' uninstalled!
	
astyle:
	@astyle --options=astyle.conf --recursive "*.h" "*.cpp"
	@./astyle-clean.sh

-include $(DEPS)
