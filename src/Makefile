all: test oalib hindices

test: bri_fortran_test.F
	f2py3 --noopt -c bri_fortran_test.F -m ftest

oalib: oalib.F
	f2py3 --noopt -c oalib.F --f77flags="-O -fPIC" --f90flags="-O -fPIC" -m oalib
	cp oalib.*.so ../oalib.so

hindices: hindices.F
	f2py3 --noopt -c hindices.F --f77flags="-O -fPIC" --f90flags="-O -fPIC" -m hindices
	cp hindices.*.so ../model/hindices.so

.PHONY: clean

clean:
	rm -rf *.so *.dSYM

