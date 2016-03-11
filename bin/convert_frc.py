#!/usr/bin/env python
"""
Simple script to convert older-style bulk forcing files with
multiple time dimensions to a new, single, unlimited time dimension
for all fields.
"""
import sys
import seapy

try:
    infile = sys.argv[1]
    outfile = sys.argv[2]
except:
    print("Usage: {:s} input_file output_file".format(sys.argv[0]))
    sys.exit()

print("Convert {:s} to {:s}".format(infile, outfile))

# Get the parameters
inc = seapy.netcdf4(infile)
lat = len(inc.dimensions['lat'])
lon = len(inc.dimensions['lon'])
epoch, tvar = seapy.roms.get_reftime(inc)

# Create the new file
onc = seapy.roms.ncgen.create_frc_bulk(
    outfile, lat=lat, lon=lon, reftime=epoch, clobber=True)

# Save the times
onc.variables['time'][:] = inc.variables[tvar][:]

# Copy the variables
fields = list(seapy.roms.forcing.fields)
fields.append('lat')
fields.append('lon')
for v in fields:
    print("{:s}, ".format(v), end='')
    onc.variables[v][:] = inc.variables[v][:]
print('done.')
onc.close()
inc.close()
