all: test oalib

test: bri_fortran_test.F
	f2py3 --noopt -c bri_fortran_test.F -m ftest

oalib: oalib.F
	f2py3 --noopt -c oalib.F --f77flags="-O -fPIC" --f90flags="-O -fPIC" -m oalib

.PHONY: clean

clean:
	rm -rf *.so *.dSYM

