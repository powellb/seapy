all: test oalib
	
test: bri_fortran_test.F
	f2py3 -c bri_fortran_test.F -m ftest

oalib: oalib.F
	f2py3 -c oalib.F --f77flags="-fPIC" -m oalib
	
