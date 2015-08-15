##State Estimation and Analysis in PYthon (SEAPY)##

Tools for working with ocean models and data.

SEAPY requires: numpy, netCDF4, joblib, and numpy_groupies

###Installation###

1) Clone the Repository into your local structure
2) Build the fortran components (you must have gfortran installed):
    % cd seapy/src && make 
3) After compilation, move the oa library into the seapy structure:
    % mv oa.so .. (or whatever your OS builds)
4) Make the documentation (you must have sphinx along with the numpydoc extension installed).
    % cd doc && make html
The documentation will then be in seapy/docs/build/html/index.html
5) If the seapy directory is not already in your python path:
    % export PYTHONPATH=`pwd`:$PYTHONPATH

You can then test by:

% python
\>\> import seapy

