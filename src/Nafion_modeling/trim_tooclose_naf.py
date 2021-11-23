#!/usr/bin/env python3
#trimming
from MolCop import mmpystream as mmps
from MolCop.analysis import topology
from MolCop.analysis import neighbor

import numpy as np
import sys
import copy

ss = mmps.Stream()
ss.import_file('dump.pos.0', 'dumppos')
print(f"Read input file ")
ss.import_file("dump.bond.0", 'dumpbond')
print(f"Read init_bond file ")

ss.sdat.create_connect_list(0.5)
connect_list = ss.sdat.connect_list
m_list = topology.create_molecule(ss.sdat)

ss.sdat.add_particles_property("molnum", _dtype=int, dim=1)
for l in m_list:
    num_mol = len(l)
    for ind in l:
        ss.sdat.particles["molnum"][ind]=num_mol
print(list(set(ss.sdat.particles["molnum"])))
mn_list_temp=list(set(ss.sdat.particles["molnum"]))
#print(ss.sdat.particles['mol'])
mn_list = [1,2,3, max(mn_list_temp)]
for i in mn_list_temp:
    if ((400 < i) & (i<430)):
        mn_list.append(i)

print("mn_list")
print(mn_list)
flag_2d = np.zeros((1,len(ss.sdat.particles['id'])), dtype=int)
for i in mn_list:
    flag = ss.sdat.particles["molnum"] == i
    print(flag)
    flag = flag.astype(int)
    flag_2d = np.append(flag_2d, [flag], axis=0)
flag1d = np.sum(flag_2d, axis=0) != 0
print(flag1d)
ss.sdat.trimming_particles(flag1d)

ss.output_file("trimmed_dump", 'dumppos', ["id","type","pos", "mask", "mol"] )

"""
nafion = copy.deepcopy(ss.sdat)
flag= (400<ss.sdat.particles['molnum']) & (ss.sdat.particles['molnum']<430)
nafion.trimming_particles(flag)
print(nafion.particles)
n_list= neighbor.make_neighbor(nafion, 1.8)

#print([i for i in n_list])
new_n_list=[]
for idx, nst in enumerate(n_list):
    #print(nst)
    for i in nst:
        if nafion.particles['mol'][idx]!=nafion.particles['mol'][i]:
            new_n_list.append([idx, nafion.particles['mol'][i]])
print(new_n_list)
"""
"""
ss2 = mmps.Stream()
ifn = sys.argv[1]
ss2.import_file(ifn, 'dumppos')
print(f"Read input file ")
ifn = sys.argv[2]
ss2.import_file(ifn, 'dumpbond')
print(f"Read init_bond file ")
ss2.sdat.create_connect_list(0.3)
connect_list = ss2.sdat.connect_list
m_list = topology.create_molecule(ss2.sdat)
ss2.sdat.add_particles_property("molnum", _dtype=int, dim=1)
for l in m_list:
    num_mol = len(l)
    for ind in l:
        ss2.sdat.particles["molnum"][ind]=num_mol
print(ss2.sdat.particles["molnum"])

flag = (ss2.sdat.particles['type']==1) & (ss2.sdat.particles["molnum"]<10)

#trimしたい短くなった炭素のidを取得する
trim_id_list=set(ss2.sdat.particles['id'][flag])
trim_id_list=list(trim_id_list)
print("trim_id_list")
print(trim_id_list)

#trimしたいidを持っている原子が入っている分子を消したい
#trimしたいidを持っている原子が入っているmolを初期構造から取り出す．
trim_mol_list = np.empty((0,1))
for i in trim_id_list:
    trim_mol= ss.sdat.particles['mol'][ss.sdat.particles['id']==i]
    print(trim_mol)
    trim_mol_list = np.append(trim_mol_list,trim_mol)
print("trim_mol_list")
print(trim_mol_list)

flag_2d = np.zeros((1,len(ss2.sdat.particles['id'])), dtype=int)
for i in trim_mol_list:
    flag = ss.sdat.particles["mol"] != i
    print(flag)
    flag = flag.astype(int) 
    flag_2d = np.append(flag_2d, [flag], axis=0)
flag1d = np.sum(flag_2d, axis=0) == len(trim_mol_list)
print(flag1d)
ss.sdat.trimming_particles(flag1d)

ss.output_file("trimmed_dump", 'dumppos', ["id","type","pos", "mask", "mol"] )
ss.output_file("trimmed_input", 'input')
"""
