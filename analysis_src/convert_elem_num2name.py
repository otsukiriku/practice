#!/usr/bin/env python3

import sys

from MolCop import mmpystream as mmps

#construct object
ss = mmps.Stream()

if len(sys.argv) <= 1:
    print("No arguments are detected")
    sys.exit()

#import input files
if len(sys.argv) >= 3:
    ifn = sys.argv[2]
    ss.import_file(ifn, 'input')
    print(f"Read input file {ifn}")

#read para.rd and set 
ss.sdat.set_elem_to_type(mmps.get_elem_to_type("para.rd"))

#import dumppos file
ifn = sys.argv[1]
#import input file
ss.import_file(ifn, 'dumppos')
print(f"Read file {ifn}")

ss.sdat.shift_particles()
ss.sdat.wrap_particles()

ss.sdat.add_particles_property("typename", _dtype=str, dim=1)
for i in ss.sdat.type_to_elem:
    flag = ss.sdat.particles['type']==i
    ss.sdat.particles['typename'][flag]=ss.sdat.type_to_elem[i]
print(ss.sdat.particles['typename'])

ofs = "ovito_" + ifn
ss.output_file(ofs, 'dumppos', ["id", "typename", "pos"])
