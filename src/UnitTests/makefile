GCC=g++
CFLAGS = -llapacke -lblas
CXXFLAGS = -g -DDEBUG
fmm: ../fmm1d_local_shift_2.cpp ../computeX_Scaled.cpp ../QR.cpp ../bsxfun.cpp testfmm.cpp
	$(GCC) $(CXXFLAGS) $^ -o $@ $(CFLAGS) 

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

testm2: ../NPart.cpp ../BinTree.cpp ../test_mat2hsssym.cpp ../QR.cpp ../compr_new.cpp ../superdcmv_cauchy.cpp ../cauchylikematvec.cpp ../superDC.cpp ../bsxfun.cpp ../superdcmv_desc.cpp ../superdcmv_node.cpp ../secular.cpp ../colnorms.cpp ../vhat.cpp ../rootfinder.cpp ../divide2.cpp testmat2hss.cpp
			$(GCC) $(CXXFLAGS) $^ -o $@ $(CFLAGS)

testmakeband: ../makeband.cpp testmakeband.cpp
			$(GCC) -DDEBUGINPUT $^ -o $@ $(CFLAGS)

testcauchylikematvec: ../cauchylikematvec.cpp ../bsxfun.cpp  test_cauchylikematvec.cpp
						g++ -DOPENBLAS -I/OpenBlas/include -std=c++11  -O3 -o cauchylikematvec.o -c ../cauchylikematvec.cpp
						g++ -DOPENBLAS -I/OpenBlas/include -std=c++11  -O3 -o test_cauchylikematvec.o -c test_cauchylikematvec.cpp
						$(GCC) -DOPENBLAS -I/OpenBlas/include -std=c++11 -L/OpenBlas/lib cauchylikematvec.o test_cauchylikematvec.o -lopenblas -o $@

test_superdcmv_cauchy: ../superdcmv_cauchy.cpp ../bsxfun.cpp ../cauchylikematvec.cpp test_superdcmv_cauchy.cpp
						g++ -DOPENBLAS -I/OpenBlas/include -std=c++11  -O3 -o cauchylikematvec.o -c ../cauchylikematvec.cpp
						g++ -DOPENBLAS -I/OpenBlas/include -std=c++11  -O3 -o test_superdcmv_cauchy.o -c test_superdcmv_cauchy.cpp
						g++ -DOPENBLAS -I/OpenBlas/include -std=c++11  -O3 -o superdcmv_cauchy.o -c ../superdcmv_cauchy.cpp
						$(GCC) -DOPENBLAS -I/OpenBlas/include -std=c++11 -L/OpenBlas/lib cauchylikematvec.o test_superdcmv_cauchy.o superdcmv_cauchy.o -lopenblas -o $@
clean:
	rm -rf fmm
	rm -rf vhat
	rm -rf colnorms
	rm -rf rootf
	rm -rf testm2
	rm -rf testmakeband
	rm -rf testcauchylikematvec
	rm -rf test_superdcmv_cauchy
	rm *.o
