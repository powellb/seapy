# !/usr/bin/env python
"""
  Module to load and parse Common Data Language (CDL) files and
  tokenize the dimensions and variables

  Written by Brian Powell on 04/30/13
  Copyright (c)2020 University of Hawaii under the MIT-License.
"""


import re


def cdl_parser(filename):
    """
    Given a netcdf-compliant CDL file, parse it to determine the structure:
    dimensions, variables, attributes, and global attributes

    Parameters
    ----------
    filename : string
        name and path of CDL file to parse

    Returns
    -------
    dims, vars, attr: dict
        dictionaries description dimensions, variables, and attributes

    """
    dim_pat = re.compile(r"\s*(\w+)\s*=\s*(\w*)\s*;")
    var_pat = re.compile(r"\s*(\w+)\s*(\w+)\({0,1}([\w\s,]*)\){0,1}\s*;")
    attr_pat = re.compile(r"\s*(\w+):(\w+)\s*=\s*\"*([^\"]*)\"*\s*;")
    global_attr_pat = re.compile(r"\s*:(\w+)\s*=\s*\"*([^\"]*)\"*\s*;")
    dims = dict()
    attr = dict()
    vars = list()
    vcount = dict()
    types = {"float": "f4", "double": "f8",
             "short": "i2", "int": "i4", "char": "S1"}

    for line in open(filename, 'r'):
        # Check if this is a dimension definition line. If it is, add
        # the dimension to the definition
        parser = dim_pat.match(line)
        if parser is not None:
            tokens = parser.groups()
            if tokens[1].upper() == "UNLIMITED":
                dims[tokens[0]] = 0
            else:
                dims[tokens[0]] = int(tokens[1])
            continue

        # Check if this is a variable definition line. If it is, add
        # the variable to the definition
        parser = var_pat.match(line)
        if parser is not None:
            tokens = parser.groups()
            nvar = {"name": tokens[1],
                    "type": types[tokens[0]],
                    "dims": tokens[2].strip().split(", ")}
            vars.append(nvar)
            vcount[tokens[1]] = len(vars) - 1
            continue

        # If this is an attribute, add the info to the appropriate variable
        parser = attr_pat.match(line)
        if parser is not None:
            tokens = parser.groups()
            if "attr" not in vars[vcount[tokens[0]]]:
                vars[vcount[tokens[0]]]["attr"] = dict()
            vars[vcount[tokens[0]]]["attr"][tokens[1]] = tokens[2]
            continue

        # If this is a global attribute, add the info
        parser = global_attr_pat.match(line)
        if parser is not None:
            tokens = parser.groups()
            attr[tokens[0]] = tokens[1]
            continue

    return dims, vars, attr


if __name__ == "__main__":
    cdl_parser("out.cdl")
