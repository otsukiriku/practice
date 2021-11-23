#!/usr/bin/env python3
import numpy as np
from MolCop import mmpystream as mmps
from MolCop.analysis import topology

ss = mmps.Stream()

print("please input 'pos file name' 'bond file name' (split with space)")

pos,bond = map(str, input().split())

ss.import_file(pos, 'dumppos')
ss.import_file(bond, 'dumpbond')

ss.sdat.set_elem_to_type(mmps.get_elem_to_type("./para.rd"))

ss.sdat.create_connect_list(0.3)
connect_list = ss.sdat.connect_list
m_list = topology.create_molecule(ss.sdat)

ss.sdat.add_particles_property("molnum", _dtype=int, dim=1)
for l in m_list:
    num_mol = len(l)
    for ind in l:
        ss.sdat.particles["molnum"][ind]=num_mol
print(ss.sdat.particles["molnum"])

n_num = ss.sdat.particles["molnum"] == 416
print("number of Nafion (which did not broken)")
print(len(ss.sdat.particles['id'][n_num])/416)

w_num = ss.sdat.particles["molnum"] == 3
print("number of water particles")
print(len(ss.sdat.particles['id'][w_num])/3)

tot = len(ss.sdat.particles['id'][n_num])+len(ss.sdat.particles['id'][w_num])

#if ss.sdat.total_particle != tot:
#    print("!! warning !! your nafionw bulk may be broken!")
