UTILITIES_DIR = ../utils
OBJ_DIR = ../../obj
BIN_DIR = ../../bin

# Includes
INCLUDES = -I$(UTILITIES_DIR)/BamTools/include/ \
           -I$(UTILITIES_DIR)/GffFileReader/ \
           -I$(UTILITIES_DIR)/Intervals/

# Source and object files
SOURCES= StatsGeneModel_main.cpp StatsGeneModel.cpp StatsGeneModel.h
OBJECTS= StatsGeneModel_main.o StatsGeneModel.o
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))

all: $(BUILT_OBJECTS)

.PHONY: all

$(BUILT_OBJECTS): $(SOURCES)
	@echo " * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/StatsGeneModel_main.o $(OBJ_DIR)/StatsGeneModel.o

.PHONY: clean
