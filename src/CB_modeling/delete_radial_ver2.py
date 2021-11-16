#!/usr/bin/env python3
import numpy as np
from MolCop import mmpystream as mmps
from MolCop.analysis import topology

ss = mmps.Stream()

ss.import_file('dump.pos.0', 'dumppos')
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

m_n = ss.sdat.particles["molnum"] == 3
print("number of water particles")
print(len(ss.sdat.particles['id'][m_n])/3)

"""
#delete H, OHterminal in inside of CB
CB_mask=[]
CB_g=[]
center=open("./center.txt", mode='r')
line = center.readline()
while True:
    line = center.readline().split()
    if not line:
        break
    CB_mask.append(line.pop(0))
    CB_g.append(line)
    print(line)
CB_mask=list(map(int, CB_mask))
#print(CB_mask)
CB_G=[]
for line in CB_g:
    CB_G.append(list(map(float, line)))
CB_G=np.array(CB_G)
#print(CB_G)
center.close()

for cent in CB_G:
    H_flag = (ss.sdat.particles['type'] != 2)|(ss.sdat.particles['type'] != 3) 
    dist = (ss.sdat.particles['pos']-cent)**2
    distarr=np.sum(dist, axis=1)
    dist_flag = distarr>9500
    flag = H_flag & dist_flag
    print(flag)
    ss.sdat.trimming_particles(flag)
    
"""

"""
for i in range(0,len(ss.sdat.particles['mol'])):
    ID = ss.sdat.particles['pos'].count_nonzero(i < 3)
print(ID)
"""
ss.output_file('withwaterlayer.dump', 'dumppos', ['id','type','pos','velo','mask'])
ss.output_file("withwaterlayer.rd", 'input')


