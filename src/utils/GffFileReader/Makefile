OBJ_DIR = ../../../obj
BIN_DIR = ../../../bin
UTILITIES_DIR = ../../utils

# Includes
INCLUDES = -I$(UTILITIES_DIR)/Intervals/

# Source and object files
SOURCES= GffFileReader.cpp GffFileReader.h
OBJECTS= $(SOURCES:.cpp=.o) 
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))

$(BUILT_OBJECTS): $(SOURCES)
	@echo " * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*

.PHONY: clean
