#!/usr/bin/env python
"""
Simple script to convert older-style bulk forcing files with
multiple time dimensions to a new, single, unlimited time dimension
for all fields.
"""
import sys
import seapy
import numpy as np


def main(infile, outfile):
    print("Convert {:s} to {:s}".format(infile, outfile))
    maxrecs = 30

    # Get the parameters
    with seapy.netcdf(infile) as inc:
        lat = len(inc.dimensions['lat'])
        lon = len(inc.dimensions['lon'])
        epoch, tvar = seapy.roms.get_reftime(inc)

        # Create the new file
        onc = seapy.roms.ncgen.create_frc_bulk(
            outfile, lat=lat, lon=lon, reftime=epoch, clobber=True)

        # Save the times
        onc.variables['time'][:] = inc.variables[tvar][:]
        ntimes = len(onc.dimensions['time'])
        onc.variables['lat'][:] = inc.variables['lat'][:]
        onc.variables['lon'][:] = inc.variables['lon'][:]

        # Copy the variables
        for v in seapy.roms.forcing.fields:
            print("{:s}, ".format(v), end='', flush=True)
            for l in seapy.chunker(np.arange(ntimes), maxrecs):
                onc.variables[v][l, :] = inc.variables[v][l, :]
        print('done.')
        onc.close()

    return 0


if __name__ == '__main__':
    try:
        infile = sys.argv[1]
        outfile = sys.argv[2]
    except:
        print("Usage: {:s} input_file output_file".format(sys.argv[0]))
        sys.exit(-1)

    sys.exit(main(infile, outfile))
