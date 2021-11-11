#!/usr/bin/env python3

from MolCop import mmpystream as mmps
from MolCop.analysis import topology

ss = mmps.Stream()

ss.import_file('show_dump', 'dumppos')
ss.import_file('dump.bond.0', 'dumpbond')

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
flag = ss.sdat.particles["molnum"] > 2
ss.sdat.trimming_particles(flag, reindex=True)

"""
for i in range(0,len(ss.sdat.particles['mol'])):
    ID = ss.sdat.particles['pos'].count_nonzero(i < 3)
print(ID)
"""
ss.output_file('test.pos', 'dumppos', ['id','type', 'mol', 'pos'])
ss.output_file("newinput", 'input')


