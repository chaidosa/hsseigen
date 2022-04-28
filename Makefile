########################################################################
####################### Makefile Template ##############################
########################################################################
OPENBLAS=0
ifeq ($(OPENBLAS), 1)
	#set BLAS_INSTALL_PATH appropriately if you have a local installation of blas library. If there is a system-wide installation available, leave this as blank
	BLAS_INSTALL_PATH=/home/nikhilh/installed_software/openblas_0_3_20
#set BLAS_LIB_NAME as blas (or mklblas or openblas or someothercustomname depending upon the library you are using).
	BLAS_LIB_NAME=openblas_nonthreaded
	CXXFLAGS = -DOPENBLAS -I$(BLAS_INSTALL_PATH)/include
	LDFLAGS =-L$(BLAS_INSTALL_PATH)/lib -l$(BLAS_LIB_NAME) 
endif


# Compiler settings - Can be customized.
CC = g++ 
CXXFLAGS += -std=c++11 
LDFLAGS = -lblas -llapacke
CXXFLAGS+= -fopenmp
#LDFLAGS+= -fopencilk
#CXXFLAGS += -DMKL_ILP64
ifeq ($(DEBUG),1)
	CXXFLAGS += -g
else
	CXXFLAGS+= -O3
endif


# Makefile settings - Can be customized.
APPNAME = Test
EXT = .cpp
SRCDIR = src
OBJDIR = obj

############## Do not change anything from here downwards! #############
SRC = $(wildcard $(SRCDIR)/*$(EXT))
OBJ = $(SRC:$(SRCDIR)/%$(EXT)=$(OBJDIR)/%.o)
#DEP = $(OBJ:$(OBJDIR)/%.o=%.d)

# UNIX-based OS variables & settings
RM = rm
DELOBJ = $(OBJ)
# Windows OS variables & settings
DEL = del
EXE = .exe
WDELOBJ = $(SRC:$(SRCDIR)/%$(EXT)=$(OBJDIR)\\%.o)
########################################################################
####################### Targets beginning here #########################
########################################################################

all: $(APPNAME)

# Builds the app
$(APPNAME): $(OBJ)
	$(CC) $(CXXFLAGS) $^ $(LDFLAGS) -o $@ 

# Creates the dependecy rules
%.d: $(SRCDIR)/%$(EXT)
	@$(CPP) $(CFLAGS) $< -MM -MT $(@:%.d=$(OBJDIR)/%.o) >$@

# Includes all .h files
#-include $(DEP) 

# Building rule for .o files and its .c/.cpp in combination with all .h
$(OBJDIR)/%.o: $(SRCDIR)/%$(EXT)
	$(CC) $(CXXFLAGS) -o $@ -c $<

################### Cleaning rules for Unix-based OS ###################
# Cleans complete project
.PHONY: clean
clean:
	$(RM) -f $(DELOBJ) $(DEP) $(APPNAME)

# Cleans only all files with the extension .d
.PHONY: cleandep
cleandep:
	$(RM) $(DEP)

#################### Cleaning rules for Windows OS #####################
# Cleans complete project
.PHONY: cleanw
cleanw:
	$(DEL) $(WDELOBJ) $(DEP) $(APPNAME)$(EXE)

# Cleans only all files with the extension .d
.PHONY: cleandepw
cleandepw:
	$(DEL) $(DEP)
