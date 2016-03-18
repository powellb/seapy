#!/usr/bin/env python
"""
Simple script to convert older-style climatology files with
multiple time dimensions to a new, single, unlimited time dimension
for all fields.
"""
import sys
import seapy
import numpy as np

try:
    infile = sys.argv[1]
    outfile = sys.argv[2]
except:
    print("Usage: {:s} input_file output_file".format(sys.argv[0]))
    sys.exit()

print("Convert {:s} to {:s}".format(infile, outfile))
maxrecs = 30

# Get the parameters
inc = seapy.netcdf(infile)
eta_rho = len(inc.dimensions['eta_rho'])
xi_rho = len(inc.dimensions['xi_rho'])
s_rho = len(inc.dimensions['s_rho'])
epoch, tvar = seapy.roms.get_reftime(inc)

# Create the new file
onc = seapy.roms.ncgen.create_clim(
    outfile, eta_rho=eta_rho, xi_rho=xi_rho, s_rho=s_rho, reftime=epoch, clobber=True)

# Save the times
onc.variables['clim_time'][:] = inc.variables[tvar][:]
ntimes = len(onc.dimensions['clim_time'])

# Copy the variables
for v in seapy.roms.fields:
    print("{:s}, ".format(v), end='', flush=True)
    for l in seapy.chunker(np.arange(ntimes), maxrecs):
        onc.variables[v][l, :] = inc.variables[v][l, :]
print('done.')
onc.close()
inc.close()
