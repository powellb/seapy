all: oalib test
	
oalib: oalib.F
	f2py -c oalib.F -m ../oalib
	
test: bri_fortran_test.F
	f2py -c bri_fortran_test.F -m ftest
