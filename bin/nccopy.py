#!/usr/bin/env python
"""
   nccopy

   Simple tool to transfer data (without attributes) from one
   file to another. It only handles variables that are in both.

   Author: Brian Powell <powellb@hawaii.edu>
   Copyright (c) 2025, Brian Powell, all rights reserved.

"""
import netCDF4
import sys
from os.path import exists


def main():
    inp = sys.argv[1]
    out = sys.argv[2]

    if not exists(inp):
        print(f"{inp} does not exist")
        sys.exit(-1)
    if not exists(out):
        print(f"{out} does not exist")
        sys.exit(-1)
        print(f"Copy data from {inp} into {out}")
    with netCDF4.Dataset(out, 'a') as of:
        with netCDF4.Dataset(inp) as inf:
            for var in inf.variables:
                if var in of.variables:
                    of.variables[var][:] = inf.variables[var][:]


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"{sys.argv[0]} srcfile destfile")
        sys.exit(-1)
    else:
        main()
