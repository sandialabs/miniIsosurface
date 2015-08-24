# -----------------------------------------------
# Mantevo Marching Cubes Utility Objects
#
# Author: sjmunn@sandia.gov
# Date  : 2015-07-31
#
# Changelog :
#   2015-07-31 - created
# -----------------------------------------------

# ========== Macros =============================

# Compiler and flags
CC              = g++
#CXXFLAGS		= -O2 -fopenmp -lprofiler -std=c++0x
#CXXFLAGS        = -O2 -lrt -std=c++0x
CXXFLAGS = -O0 -g -lrt -std=c++0x

# Source and build directories
SRCDIRS         = Data_Obj IO Reporting Constants User_Interface Algorithm GeneralContext
BUILD_DIRECTORY = build
BUILDDIR        = $(BUILD_DIRECTORY)

# Get source code files and generate object names
SOURCES        := $(foreach dir,$(SRCDIRS),$(wildcard $(dir)/*.cpp))
OBJECTS        := $(addprefix $(BUILDDIR)/,$(notdir $(SOURCES:.cpp=.o)))

# ========== Compiler Rules =====================

# Linker commands

all: $(OBJECTS)
	@echo "All utility objects have been compiled";
	
# General purpose object compiler
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
