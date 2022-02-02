GCC=g++
CFLAGS = -llapacke -lcblas -llapack
#fmm: ../fmm1d_local_shift.cpp testfmm.cpp
#	$(GCC) $(CFLAGS) $^ -o $@

#vhat: ../vhat.cpp testvhat.cpp
#	$(GCC) $(CFLAGS) $^ -o $@

#colnorms: ../bsxfun.cpp ../colnorms.cpp test_colnrm.cpp
#		$(GCC)  $^ -o $@ $(CFLAGS)

rootf: ../bsxfun.cpp ../rootfinder.cpp testRootfinder.cpp
		$(GCC)  $^ -o $@ $(CFLAGS)

clean:
	rm -rf fmm
	rm -rf vhat
	rm -rf colnorms
	rm -rf rootf