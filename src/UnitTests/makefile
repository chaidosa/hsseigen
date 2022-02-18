GCC=g++
CFLAGS = -llapacke -lcblas -llapack
CXXFLAGS = -g
fmm: ../fmm1d_local_shift.cpp testfmm.cpp
	$(GCC) $(CFLAGS) $^ -o $@

vhat: ../vhat.cpp testvhat.cpp
	$(GCC) $(CFLAGS) $^ -o $@

colnorms: ../bsxfun.cpp ../colnorms.cpp test_colnrm.cpp
		$(GCC)  $^ -o $@ $(CFLAGS)

rootf: ../bsxfun.cpp ../rootfinder.cpp testRootfinder.cpp
		$(GCC)  $^ -o $@ $(CFLAGS)

tb2hss: ../BinTree.cpp ../NPart.cpp ../band2hss.cpp testBand2hss.cpp
		$(GCC) $(CXXFLAGS) $^ -o $@ $(CFLAGS)

testsec: ../rootfinder.cpp ../bsxfun.cpp ../secular.cpp ../colnorms.cpp ../vhat.cpp testSEC.cpp
		$(GCC) $(CXXFLAGS) $^ -o $@ $(CFLAGS)

testsuperdcmv: ../superdcmv_cauchy.cpp ../cauchylikematvec.cpp ../bsxfun.cpp test_superdcmvs.cpp 
				$(GCC) $(CXXFLAGS) $^ -o $@ $(CFLAGS)

testm2: ../NPart.cpp ../BinTree.cpp ../test_mat2hsssym.cpp ../QR.cpp ../compr_new.cpp ../superdcmv_cauchy.cpp ../cauchylikematvec.cpp ../superDC.cpp ../bsxfun.cpp ../superdcmv_desc.cpp ../superdcmv_node.cpp ../secular.cpp ../colnorms.cpp ../vhat.cpp ../rootfinder.cpp ../divide2.cpp teshmat2hss.cpp
			$(GCC) $(CXXFLAGS) $^ -o $@ $(CFLAGS)
clean:
	rm -rf fmm
	rm -rf vhat
	rm -rf colnorms
	rm -rf rootf
	rm -rf testm2