########################################################################
####################### Makefile Template ##############################
########################################################################
OPENBLAS=1
ifeq ($(OPENBLAS), 1)
	#set BLAS_INSTALL_PATH appropriately if you have a local installation of blas library. If there is a system-wide installation available, leave this as blank
	# BLAS_INSTALL_PATH=/OpenBlas
#set BLAS_LIB_NAME as blas (or mklblas or openblas or someothercustomname depending upon the library you are using).
	# BLAS_LIB_NAME=openblas
	# CXXFLAGS = -DOPENBLAS -I$(BLAS_INSTALL_PATH)/include
	CXXFLAGS = -DOPENBLAS
	LDFLAGS =-lblas -llapack -llapacke 
	# LDFLAGS = -L$(BLAS_INSTALL_PATH)/lib -l$(BLAS_LIB_NAME)
	CC=g++
else
# 	CC=icpx
# 	LDFLAGS = -lmkccl_rt
	LDFLAGS =-lblas -llapack -llapacke 
	CC=g++
endif

# Compiler settings - Can be customized.
CXXFLAGS += -std=c++11 
ifeq ($(PARALLEL),1)
	CXXFLAGS += -fopenmp -DPARALLEL
endif
#LDFLAGS+= -fopencilk
#CXXFLAGS += -DMKL_ILP64
ifeq ($(DEBUG),1)
	CXXFLAGS += -g -DDEBUG
else
	CXXFLAGS += -O3
endif

# Makefile settings - Can be customized.
APPNAME = Test
EXT = .cpp
SRCDIR = src
OBJDIR = obj

SRC = $(wildcard $(SRCDIR)/*$(EXT))
OBJ = $(SRC:$(SRCDIR)/%$(EXT)=$(OBJDIR)/%.o)

all: $(APPNAME)

# Builds the app
$(APPNAME): $(OBJ)
	$(CC) $(CXXFLAGS) $^ $(LDFLAGS) -o Test	

# Creates the dependecy rules
#%.d: $(SRCDIR)/%$(EXT)
#	@$(CPP) $(CFLAGS) $< -MM -MT $(@:%.d=$(OBJDIR)/%.o) >$@

# Includes all .h files
#-include $(DEP) 

# Building rule for .o files and its .c/.cpp in combination with all .h
$(OBJDIR)/%.o: $(SRCDIR)/%$(EXT)
	$(CC) $(CXXFLAGS) -o $@ -c $<

# Cleans complete project
.PHONY: clean
clean:
	rm -rf $(OBJ) $(DEP) $(APPNAME)

# Cleans only all files with the extension .d
.PHONY: cleandep
cleandep:
	rm $(DEP)
