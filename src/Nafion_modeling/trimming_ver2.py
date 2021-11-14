#!/usr/bin/env python3
#trimming
from MolCop import mmpystream as mmps
from MolCop.analysis import topology

import numpy as np
import sys
import copy

CB_num = 4
Pt_num = 56
def g_c_by_mask(atoms : mmps.Stream().sdat.particles, mask):
    flag = atoms["mask"] == mask
    center = atoms["pos"][flag].mean(axis=0)
    return center

ss = mmps.Stream()
ifn = sys.argv[1]
ss.import_file(ifn, 'input')
print(f"Read input file {ifn}")
ss.import_file("dump.pos.20000", 'dumppos')
print(f"Read dump file ")

ss.sdat.wrap_particles()
ss.sdat.shift_particles()

ss.import_file("dump.bond.20000", 'dumpbond')
print(f"Read bond file ")

ss.sdat.create_connect_list(0.3)

connect_list = ss.sdat.connect_list
m_list = topology.create_molecule(ss.sdat)

flag = ss.sdat.particles["mol"] == 1
ss.sdat.trimming_particles(flag, reindex=True)

#read_CBcenter
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

Type=ss.sdat.particles['type']
#delete H, OHterminal in inside of CB
for cent in CB_G:
    H_flag = ((Type != 2)|(Type != 3)) 
    dist = (ss.sdat.particles['pos']-cent)**2
    distarr=np.sum(dist, axis=1)
    dist_flag = distarr>10000
    flag = H_flag | dist_flag
    print(flag)
    ss.sdat.trimming_particles(flag)
    
Pt_list=[i for i in range(7,7+Pt_num)]
Pt_G = np.array([g_c_by_mask(ss.sdat.particles, l) for l in Pt_list])

elem_num=len(ss.sdat.particles['id'])
flag_2d_CB = np.zeros((1,elem_num), dtype=int)
for c in CB_G:
    dist = (ss.sdat.particles['pos'] - c)**2
    #同じ行のものを足す
    distarr=np.sum(dist, axis=1)
    flag = distarr<11000
    #bool -> int 
    flag = flag.astype(int)
    flag_2d_CB = np.append(flag_2d_CB, [flag], axis=0)
dist_flag1 = np.sum(flag_2d_CB, axis=0) > 0


flag_2d_Pt = np.zeros((1,elem_num), dtype=int)
for pt_g in Pt_G:
    dist = (ss.sdat.particles['pos']-pt_g)**2
    distarr=np.sum(dist, axis=1)
    dist_flag2 = distarr > 14*14
    flag_2d_Pt = np.append(flag_2d_Pt, [dist_flag2], axis=0)
dist_flag2 = np.sum(flag_2d_Pt.astype(int), axis=0)
print(dist_flag2)
print(dist_flag2.max())
dist_flag2 = (dist_flag2 == dist_flag2.max()) | (Type==4)  

#flag = dist_flag2
#print(flag)
flag = dist_flag1 | dist_flag2
    
ss.sdat.trimming_particles(flag)


ss.output_file("trimmed_dump", 'dumppos', ["id","type","pos","mask"] )
ss.output_file("trimmed_input", 'input')

