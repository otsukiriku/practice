#!/usr/bin/env python3
import numpy as np
import copy
from MolCop import mmpystream as mmps
from MolCop.analysis import topology

ss = mmps.Stream()
ss_out = mmps.Stream()

ss.import_file('dump.pos.0', 'dumppos')
ss.import_file('dump.bond.0', 'dumpbond')

ss.sdat.set_elem_to_type(mmps.get_elem_to_type("./para.rd"))

"""
ss.sdat.create_connect_list(0.3)
connect_list = ss.sdat.connect_list
m_list = topology.create_molecule(ss.sdat)

ss.sdat.add_particles_property("molnum", _dtype=int, dim=1)
for l in m_list:
    num_mol = len(l)
    for ind in l:
        ss.sdat.particles["molnum"][ind]=num_mol

print(len(ss.sdat.particles['id'][ss.sdat.particles["molnum"]==415]))

water = ss.sdat.particles["mol"]
ss_out.sdat.cell = ss.sdat.cell
ss_out.sdat.dcell = ss.sdat.dcell
ss_out.sdat.newcell = ss.sdat.newcell
water = ss.sdat.particles
ss_out.sdat.concate_particles(ss.sdat.particles)
ss_out.sdat.concate_particles({'type':ss.sdat.particles['type'][flag]})
ss_out.sdat.concate_particles({'pos':ss.sdat.particles['pos'][flag]})
"""

ss.sdat.sort_particles_by_id()
ss_out.sdat.cell = ss.sdat.cell
ss_out.sdat.dcell = ss.sdat.dcell
ss_out.sdat.newcell = ss.sdat.newcell
flag=ss.sdat.particles['id'] > 89856
print(ss.sdat.particles['id'][flag])
water = copy.deepcopy(ss.sdat)
water.trimming_particles(flag)
ss_out.sdat.concate_particles(water.particles)

flag=ss.sdat.particles['id'] < 89857
ss.sdat.trimming_particles(flag)
ss_out.sdat.concate_particles(ss.sdat.particles)

"""
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

m_n = ss.sdat.particles["molnum"] == 3
print("number of water particles")
print(len(ss.sdat.particles['id'][m_n])/3)
"""

ss_out.output_file('newnafionw.dump', 'dumppos', ['id','type','pos','velo'])
#ss_out.output_file("withwaterlayer.rd", 'input')


