UTILITIES_DIR = ../utils
OBJ_DIR = ../../obj
BIN_DIR = ../../bin

# Includes
INCLUDES = -I$(UTILITIES_DIR)/BamTools/include/ \
	   -I$(UTILITIES_DIR)/GffFileReader/ \
           -I$(UTILITIES_DIR)/Intervals/

# Source and object files
SOURCES= CountGeneModel_main.cpp CountGeneModel.cpp CountGeneModel.h 
OBJECTS= CountGeneModel_main.o CountGeneModel.o 
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))

all: $(BUILT_OBJECTS)

.PHONY: all

$(BUILT_OBJECTS): $(SOURCES)
	@echo " * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/CountGeneModel_main.o $(OBJ_DIR)/CountGeneModel.o

.PHONY: clean
