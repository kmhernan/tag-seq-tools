OBJ_DIR = ../../obj
BIN_DIR = ../../bin

# Source and object files
SOURCES= CountMatrix_main.cpp CountMatrix.cpp CountMatrix.h
OBJECTS= CountMatrix_main.o CountMatrix.o
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))

all: $(BUILT_OBJECTS)

.PHONY: all

$(BUILT_OBJECTS): $(SOURCES)
	@echo " * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/CountMatrix_main.o $(OBJ_DIR)/CountMatrix.o

.PHONY: clean
