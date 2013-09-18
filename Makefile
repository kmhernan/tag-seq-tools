# =============================================
# TagSeqTools Makefile
# Kyle Hernandez, UNLICENSED 2013
# The University of Texas at Austin
# Juenger Lab
# kmhernan84@gmail.com
#
# adapted from bedtools, (c) 2009 Aaron Quinlan
# =============================================

SHELL := /bin/bash -e

# check whether the user is on Linux or OS X
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
	export CXX = clang++ 
	export CXXFLAGS +=-std=c++11 -stdlib=libc++ -Wall -O2 -D_FILE_OFFSET_BITS=64 -fPIC
endif
ifeq ($(UNAME_S),Linux)
	export CXX = g++
	export CXXFLAGS =-std=gnu++11 -Wall -O2 -D_FILE_OFFSET_BITS=64 -fPIC
endif

# define our object and binary directories
export OBJ_DIR	= obj
export BIN_DIR	= bin
export SRC_DIR	= src
export UTIL_DIR	= src/utils
export LIBS     = -lz
export BT_ROOT  = src/utils/BamTools/

SUBDIRS = $(SRC_DIR)/CountGeneModel \
	  $(SRC_DIR)/CountMatrix

UTIL_SUBDIRS =	$(SRC_DIR)/utils/BamTools \
	     	$(SRC_DIR)/utils/gzstream \
                $(SRC_DIR)/utils/GffFileReader

BUILT_OBJECTS = $(OBJ_DIR)/*.o

all: print_banner $(OBJ_DIR) $(BIN_DIR) $(UTIL_SUBDIRS) $(SUBDIRS)
	@echo "- Building main TagSeqTools binary."
	@$(CXX) $(CXXFLAGS) -c src/TagSeqTools.cpp -o obj/TagSeqTools.o
	@$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $(BIN_DIR)/TagSeqTools $(BUILT_OBJECTS) -L$(UTIL_DIR)/BamTools/lib/ -lbamtools $(LIBS)
	@echo "done."

.PHONY: all

print_banner:
	@echo "Building TagSeqTools:"
	@echo "================================================"
.PHONY: print_banner

# make the "obj/" and "bin/" directories, if they don't exist
$(OBJ_DIR) $(BIN_DIR):
	@mkdir -p $@


# One special case: All (or almost all) programs requires the BamTools API files to be created first.
.PHONY: bamtools_api
bamtools_api:
	@$(MAKE) --no-print-directory --directory=$(BT_ROOT) api
$(UTIL_SUBDIRS) $(SUBDIRS): bamtools_api


# even though these are real directories, treat them as phony targets, forcing to always go in them are re-make.
# a future improvement would be the check for the compiled object, and rebuild only if the source code is newer.
.PHONY: $(UTIL_SUBDIRS) $(SUBDIRS)
$(UTIL_SUBDIRS) $(SUBDIRS): $(OBJ_DIR) $(BIN_DIR)
	@echo "- Building in $@"
	@$(MAKE) --no-print-directory --directory=$@

clean:
	@$(MAKE) --no-print-directory --directory=$(BT_ROOT) clean_api
	@echo " * Cleaning up."
	@rm -rf $(OBJ_DIR)/ $(BIN_DIR)/
.PHONY: clean
