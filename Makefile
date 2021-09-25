CPP=g++ #-std=c++11
BLASHOME = ./OpenBLAS
#INCLUDELIB = $(BLASHOME)/bin
INCLUDELIB = $(BLASHOME)/lib
INCLUDES += -I$(BLASHOME)/include -I.
PAPI_HOME=/home/min/a/hegden/Research/PAPI
ifeq ($(DEBUG),1)
CFLAGS=-g -DDEBUG
else
CFLAGS=-DNDEBUG
endif

ifeq ($(PAPI),1)
CFLAGS += -DPAPI
CFLAGS +=$(PAPI_HOME)/lib/libpapi.so
endif

BIN_DEPS = Test.cpp QR.cpp mat2hsssym.cpp BinTree.cpp makeband.cpp NPart.cpp compr.cpp RandGen.cpp Divide.cpp

#INCLUDES += $(INCLUDELIB)/libopenblas.dll
INCLUDES += $(INCLUDELIB)/libopenblas_haswellp-r0.3.10.so#libopenblas.so

all : hsseigen

hsseigen: $(BIN_DEPS)
	$(CPP) $(CFLAGS) -o $@ $^ $(INCLUDES) 
clean:
	-rm hsseigen.exe hsseigen
