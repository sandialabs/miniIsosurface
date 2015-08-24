# -----------------------------------------------
# Mantevo Marching Cubes OpenMP Implementation
# with Merger algorithm
#
# Author: sjmunn@sandia.gov
# Date  : 2015-07-31
#
# Changelog :
#   2015-07-31 - created
# -----------------------------------------------

# ========== Macros =============================

# project name (generate executable with this name)
TARGET          = MCMpi.exe

# Compiler and flags
CC              = mpic++
#CXXFLAGS        = -O2 -fopenmp -std=c++0x
CXXFLAGS        = -O0 -g -fopenmp -std=c++0x

# Source and build directories
SRCDIRS         = Algorithm Implementations User_Interface
BUILD_DIRECTORY = build
BUILDDIR        = $(BUILD_DIRECTORY)

# Utility object code
UTILS_DIR       = ../common/build

# Get source code files and generate object names
SOURCES        := $(foreach dir,$(SRCDIRS),$(wildcard $(dir)/*.cpp))
OBJECTS        := $(addprefix $(BUILDDIR)/,$(notdir $(SOURCES:.cpp=.o)))
UTILS_OBJ	   := $(wildcard $(UTILS_DIR)/*.o)

# ========== Compiler Rules =====================

# Linker commands

all: main.o $(OBJECTS) $(UTILS_OBJ)
	@echo "Linking object files...";
	@$(CC) $(CXXFLAGS) $(BUILDDIR)/main.o $(OBJECTS) $(UTILS_OBJ) -o $(BUILDDIR)/$(TARGET) $(GPERF_FLAG);
	@echo "Compiled and linked "$(BUILDDIR)/$(TARGET)" successfully!"

# Compiler commands

main.o: main.cpp
	@mkdir -p $(BUILDDIR);
	@echo "Compiling "$(BUILDDIR)/$@;
	@$(CC) $(CXXFLAGS) -c main.cpp -o $(BUILDDIR)/$@;

# General purpose object code compiler
$(BUILDDIR)/%.o: ./*/%.cpp
	@mkdir -p $(BUILDDIR);
	@echo "Compiling "$@;
	@$(CC) $(CXXFLAGS) -c $< -o $@;

# ========== Cleanup Rules ======================

.PHONEY: clean
clean:
	rm -f $(BUILDDIR)/*.o;
	
.PHONEY: clean_run_data
clean_run_data:
	rm -f yaml_out.*;
	rm -f callgrind.*;

.PHONEY: deepclean
deepclean: clean clean_run_data
	rm -rf $(BUILDDIR);
	
# ========== Debug Rules ========================

debug:
	@echo "Using compiler: "
	@echo $(CC)
	@echo "With flags: "
	@echo $(CXXFLAGS)
	@echo "Source code directories: "
	@echo $(SRCDIRS);
	@echo "Source code files: "
	@echo $(SOURCES);
	@echo "Object files: "
	@echo $(OBJECTS);
