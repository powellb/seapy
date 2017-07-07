#!/usr/bin/env python

import seapy
import ipdb

# Call the code we want to debug here... Make sure to insert:
# ipdb.set_trace()
# where we wish to begin our debugging trace
ipdb.set_trace()
t = seapy.roms.obsgen.modis_sst_map("grid.nc", 5)
