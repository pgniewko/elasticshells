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
SOURCES	     := main.cpp \
		$(shell find $(SRC) -type f -name "*.cpp")

HEADERS	     := $(shell find $(SRC) -type f -name "*.h")

OBJECTS      := $(SOURCES:.cpp=.o)

TEST_OBJECTS := $(TEST_SOURCES:.cpp=.o)

DEPS         := $(OBJECTS:.o=.d)

#Linking commands:
$(TARGET): $(OBJECTS)
	@echo LINKING ...
	$(MKDIR) $(@D)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)
	@echo BUILDING IS DONE

#Compilation commands:
main.o: $(SOURCES) $(HEADERS) #main.cpp 

$(SRC)/%.o: $(SRC)/%.cpp $(SRC)/%.h

%.o: %.cpp %.h
	$(CXX) -MMD -MP $(CXXFLAGS) -c $< -o $@
	
# Tell make that these are phony targets
.PHONY: build clean tests install uninstall indent


build: $(TARGET)
	@echo BUILDING IS DONE

clean:
	@echo Cleaning...
	$(RM) $(TARGET) $(OBJECTS) $(DEPS)
	
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
