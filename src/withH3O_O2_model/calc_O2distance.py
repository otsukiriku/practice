#!/usr/bin/env python3
from MolCop import mmpystream as mmps
from MolCop.analysis import topology

import sys
import numpy as np
import math



ss = mmps.Stream()
ss.import_file(sys.argv[1], 'dumppos')
ss.import_file(sys.argv[2], 'dumpbond')
print(f"Read dump file ")
ss.sdat.sort_particles_by_id();
ss.sdat.create_connect_list();
#print(ss.sdat.particles)
Pt_id = ss.sdat.particles["id"][ss.sdat.particles["type"]==4]
#Pt_end_id = Pt_id[-1]

O2_id = ss.sdat.particles["id"][ss.sdat.particles["type"]==16]
#print(ss.sdat.connect_list[111145])
O2_pos = ss.sdat.particles["pos"][ss.sdat.particles["type"]==16]
O2_start_id = O2_id[0]
#O2_connect_list = 
O2_visited_flag = [False]*(len(O2_id))
O_connect_with_Pt =  [False]*(len(O2_id))

O2_dist_connect_withO = []
O2_bondorder_connect_withO = []
O2_dist_NOT_connect_withO = []
O2_bondorder_NOT_connect_withO = []


for current_id,Oid in enumerate(O2_id):
    if O2_visited_flag[current_id]:
        continue

    O2_visited_flag[current_id] = True
    tmp_connect_list = ss.sdat.connect_list[Oid]
    tmp_bondorder_list = ss.sdat.bondorder_list[Oid]
    for t_c_l,t_b_l in zip(tmp_connect_list,tmp_bondorder_list):
        pair_id = t_c_l-O2_start_id
        if t_c_l in Pt_id:
            O_connect_with_Pt[current_id] = True
        if (t_c_l in O2_id) and O2_visited_flag[pair_id]:
            tmp_dist = np.sum( (ss.sdat.particles["pos"][Oid]-ss.sdat.particles["pos"][t_c_l])**2, axis=0)
            tmp_dist = math.sqrt(tmp_dist.tolist())
            if O_connect_with_Pt[current_id] or O_connect_with_Pt[pair_id]:
                O2_dist_connect_withO.append(tmp_dist)
                O2_bondorder_connect_withO.append(t_b_l[-1])
            else:
                if tmp_dist<2:
                    O2_dist_NOT_connect_withO.append(tmp_dist)
                    O2_bondorder_NOT_connect_withO.append(t_b_l[-1])

print("O-O of O2 molecule distance on Pt = {}".format(sum(O2_dist_connect_withO)/len(O2_dist_connect_withO)))
print("O-O of O2 molecule bondorder on Pt = {}".format(sum(O2_bondorder_connect_withO)/len(O2_bondorder_connect_withO)))
print("O-O of O2 molecule distance NOT on Pt = {}".format(sum(O2_dist_NOT_connect_withO)/len(O2_dist_NOT_connect_withO)))
print("O-O of O2 molecule distance NOT on Pt = {}".format(sum(O2_bondorder_NOT_connect_withO)/len(O2_bondorder_NOT_connect_withO)))



