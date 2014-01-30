all: test oalib
	
test: bri_fortran_test.F
	f2py -c bri_fortran_test.F -m ftest

oalib: oalib.F
	f2py -c oalib.F --f77flags="-fPIC" -m oalib
	
