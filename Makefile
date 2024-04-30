########################################################################
####################### Makefile Template ##############################
########################################################################
ifeq ($(OPENBLAS), 1)
	#set BLAS_INSTALL_PATH appropriately if you have a local installation of blas library. If there is a system-wide installation available, leave this as blank
	BLAS_INSTALL_PATH=/home/c200010021/software/openblasinstall
#set BLAS_LIB_NAME as blas (or mklblas or openblas or someothercustomname depending upon the library you are using).
	BLAS_LIB_NAME=openblas
	CXXFLAGS = -DOPENBLAS -I$(BLAS_INSTALL_PATH)/include
	# CXXFLAGS = -DOPENBLAS
	LDFLAGS =-lblas -llapack -llapacke 
	LDFLAGS = -L$(BLAS_INSTALL_PATH)/lib -l$(BLAS_LIB_NAME)
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

ifeq ($(TASK),1)
	CXXFLAGS += -fopenmp -DPARALLEL_TASK
endif

ifeq ($(CILK),1)
	CC = clang++
	CXXFLAGS += -fopencilk -DOPEN_CILK
endif
#LDFLAGS+= -fopencilk
#CXXFLAGS += -DMKL_ILP64
ifeq ($(DEBUG),1)
	CXXFLAGS += -g -DDEBUG
else
	CXXFLAGS += -O3
endif

ifeq ($(DIST), 1)
	CC = mpic++
	CXXFLAGS += -DDIST
endif

ifeq ($(HYBRD), 1)
	CC = mpic++
	CXXFLAGS += -DHYBRD -fopenmp
endif

# Makefile settings - Can be customized.
APPNAME = Test
EXT = .cpp
SRCDIR = src
OBJDIR = obj

SRC = $(wildcard $(SRCDIR)/*$(EXT))
OBJ = $(SRC:$(SRCDIR)/%$(EXT)=$(OBJDIR)/%.o)

all: objdir $(APPNAME)

# Builds the app
$(APPNAME): $(OBJ)
	$(CC) $(CXXFLAGS) $^ $(LDFLAGS) -o Test	

# Creates the dependecy rules
#%.d: $(SRCDIR)/%$(EXT)
#	@$(CPP) $(CFLAGS) $< -MM -MT $(@:%.d=$(OBJDIR)/%.o) >$@

# Includes all .h files
#-include $(DEP) 

.PHONY: objdir
objdir:
	mkdir -p obj/

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

CILKBENCH = Testbench
CILKSCALE = Testscale

objbench:
	rm -rf obj/bench/
	mkdir -p obj/bench/

objscale:
	rm -rf obj/scale/
	mkdir -p obj/scale/

BENCHAPP = Testbench

OBJBench = $(SRC:$(SRCDIR)/%$(EXT)=$(OBJDIR)/bench/%.o)

$(OBJDIR)/bench/%.o: $(SRCDIR)/%$(EXT)
	clang++ $(CXXFLAGS) -fopencilk -DOPEN_CILK -fcilktool=cilkscale-benchmark -o $@ -c $<

$(BENCHAPP): $(OBJBench)
	clang++ $(CXXFLAGS) -fopencilk -DOPEN_CILK -fcilktool=cilkscale-benchmark $^ $(LDFLAGS) -o Testbench

Benchmark: objbench $(BENCHAPP)


SCALEAPP = Testscale

OBJScale = $(SRC:$(SRCDIR)/%$(EXT)=$(OBJDIR)/scale/%.o)

$(OBJDIR)/scale/%.o: $(SRCDIR)/%$(EXT)
	clang++ $(CXXFLAGS) -fopencilk -DOPEN_CILK -fcilktool=cilkscale -o $@ -c $<

$(SCALEAPP): $(OBJScale)
	clang++ $(CXXFLAGS) -fopencilk -DOPEN_CILK -fcilktool=cilkscale $^ $(LDFLAGS) -o Testscale

Scale: objscale $(SCALEAPP)

runSortScale: # Scale Benchmark
	python3 /home/sysad/customsoftware/opencilk_2_0_0/share/Cilkscale_vis/cilkscale.py -c ./Testscale -b ./Testbench -cpus 1,2,3,4,8,16,32,48  --output-csv report.csv --output-plot report.pdf --args ../hssdata/triad_8k.txt 65536 64 2 1 96 10
