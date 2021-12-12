#!/usr/bin/env python3
from MolCop import mmpystream as mmps
from MolCop.analysis import topology
from E_T.func import read_center as rc
import numpy as np
import sys
import copy
sys.path.append("/nfshome12/rotsuki/molcop/src/")
#print(sys.path)
import modeling.unwrap_particles as u_p
import evaluate_structure.get_center as g_c

Pt_cluster_num=326
Pt_num=14
#cut off distance from Pt surface
cutoff = 9
target = 5
# 5:SofNafion 6:FofNadion

CB_center = rc.r_c("new_center.txt")
print(CB_center)

ss = mmps.Stream()
ifn = sys.argv[1]
ss.import_file(ifn, 'dumppos')
ifn = sys.argv[2]
ss.import_file(ifn, 'dumpbond')
print(f"Read dump file ")

#Get Ptcenter
#Pt_G = gpc.get_Pt_G(ss.sdat, Pt_num, Pt_cluster_num)

ss.sdat.create_connect_list(0.3)
connect_list = ss.sdat.connect_list

#flagconnect reindexの両方をTrueにすることで
#connect_listを残してトリミングすることもできる．
Pt = copy.deepcopy(ss.sdat)
flag = ss.sdat.particles['type']==4
Pt.flagconnect = True
Pt.trimming_particles(flag, reindex=True)

#配位数でPtの表面を管理することにした．
coordination_num=[]
coordination_num.extend([len(i) for i in Pt.connect_list])
coordination_num = np.array(coordination_num)

distflag = np.zeros_like(Pt.particles['id'], dtype=np.bool)
#print(distflag)
for cent in CB_center:
    dist = (Pt.particles['pos'] - cent)**2
    dist = np.sum(dist, axis=1)
    #print(*dist)
    distarr = dist < 10800 
    distflag += distarr
#flag反転し，中心から一定以上離れたものを残す．
distflag = np.logical_not(distflag)
flag = (coordination_num < 10) & distflag
Pt.trimming_particles(flag, reindex=True)
m_list = topology.create_molecule(Pt)
Pt_mol=list(set(Pt.particles['mol']))
#Ptはこれで表面だけ残った．
#trim only target particle
taeget_flag = ss.sdat.particles['type'] == target
#F_flag = ss.sdat.particles['type'] == 6
ss.sdat.trimming_particles(taeget_flag)

target = ss.sdat

count=[]
for p_mol in Pt_mol:
    pt_pos = Pt.particles['pos'][Pt.particles['mol'] == p_mol]
    distflag = np.zeros_like(target.particles['id'], dtype=np.bool)
    for p in pt_pos:
        dist = (target.particles['pos'] - p)**2
        dist = np.sum(dist, axis=1)
        flag = dist < (cutoff**2)
        distflag += flag
    temp=np.sum(distflag)
    temp.tolist()
    count.append(temp)

print('Pt_mol_ID Number_of_target_atom')
for idx, c in enumerate(count):
    print(idx+1, c)


#for pt_pos in 

#ss1 = mmps.Stream()
#ss1.sdat = Pt

#ss1.output_file("show_dump", 'dumppos', ["id","type","pos", "mol"] )
#ss.output_file("trimmed_input", 'input')

