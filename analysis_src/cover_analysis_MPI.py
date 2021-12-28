#!/usr/bin/env python3
from MolCop import mmpystream as mmps
from MolCop.analysis import topology
from E_T.func import read_center as rc
from E_T.func import get_center as gc
import numpy as np
import sys
import copy
import math
sys.path.append("/nfshome12/rotsuki/molcop/src/")
#print(sys.path)
import modeling.unwrap_particles as u_p
import evaluate_structure.get_center as g_c

#for MPI
import os
from mpi4py import MPI


#********************************
Pt_cluster_num=309
Pt_num=56
Pt_mask_start = 7
#cut off distance from Pt surface
cutoff = 9
target_type = 6
# 5:SofNafion 6:FofNadion
#********************************



def distance_flag(atoms: mmps.Stream().sdat, center=[]):
    distflag = np.zeros_like(atoms.particles['id'], dtype=bool)
    #print(distflag)
    for cent in center:
        dist = (atoms.particles['pos'] - cent)**2
        dist = np.sum(dist, axis=1)
        #print(*dist)
        distarr = dist < 10800
        distflag += distarr
    #flag反転し，中心から一定以上離れたものを残す．
    distflag = np.logical_not(distflag)
    return distflag

def get_coordination(atoms: mmps.Stream().sdat):
    coordination_num=[]
    coordination_num.extend([len(i) for i in atoms.connect_list])
    coordination_num = np.array(coordination_num)
    return coordination_num

def calc_coverage(atoms:mmps.Stream().sdat, target:mmps.Stream().sdat, conditions=[]):
    count=[]
    cover=[0 for _ in conditions]
    for idx,condtion in enumerate(conditions):
        pt_pos = atoms.particles['pos'][atoms.particles['mask'] == condtion]
        distflag = np.zeros_like(target.particles['id'], dtype=bool)
        for p in pt_pos:
            dist = (target.particles['pos'] - p)**2
            dist = np.sum(dist, axis=1)
            #[print(math.sqrt(d)) for d in dist if d < cutoff**2]
            flag = dist < (cutoff**2)
            sumflag = np.sum(flag)
            coverflag=0
            if sumflag != 0:
                coverflag = 1
            cover[idx]+=coverflag
            distflag += flag
        temp=np.sum(distflag)
        temp.tolist()
        count.append(temp)
    return count,cover

def main():
    # for MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    if rank == 0:
        print(f"Input process : {size}")
    CB_center = rc.r_c("./center.txt")
    #print(CB_center)
    ss = mmps.Stream()
    ss.import_file('config.rd', 'config')

    # MPI process
    start_step = int(sys.argv[1])
    end_step = int(sys.argv[2])
    total_step = end_step - start_step

    file_step = ss.sdat.file_step
    dstep = int(total_step/(file_step*size)) + 1

    init_step = dstep * rank * file_step + start_step
    final_step = min(dstep * (rank+1) * file_step + start_step, total_step + start_step)

    if final_step == end_step:
        final_step += 1


    sumcount = []
    sumcover = []
    for step in range(init_step, final_step, file_step):
        ifn = 'dump.pos.' + str(step)
        if not os.path.isfile(ifn):
            print(f"{ifn} is not exit")
            continue
        ss.import_file(ifn, 'dumppos')
        ifn = 'dump.bond.' + str(step)
        if not os.path.isfile(ifn):
            print(f"{ifn} is not exit")
            continue
        ss.import_file(ifn, 'dumpbond')
        #Get surface Pt info
        ss.sdat.create_connect_list(0.3)
        connect_list = ss.sdat.connect_list
        #flagconnect reindexの両方をTrueにすることで
        #connect_listを残してトリミングすることもできる．
        Pt = copy.deepcopy(ss.sdat)
        flag = ss.sdat.particles['type']==4
        Pt.flagconnect = True
        Pt.trimming_particles(flag, reindex=True)

        #配位数でPtの表面を管理することにした．
        coordination_num = get_coordination(Pt)
        distflag = distance_flag(Pt, CB_center)
        flag = (coordination_num < 10) & distflag
        Pt.trimming_particles(flag, reindex=True)
        topology.unwrap_molecule(Pt)
        #m_list = topology.create_molecule(Pt)
        Pt.wrap_particles()
        #Ptはこれで表面だけ残った．
        Pt_mask = [i for i in range(Pt_mask_start, Pt_mask_start+Pt_num)]

        #print("Now processing dump.pos."+str(ifn))
        #trim only target particle
        target_flag = ss.sdat.particles['type'] == target_type
        ss.sdat.trimming_particles(target_flag)
        Target = ss.sdat

        #print(cover)
        count, cover = calc_coverage(Pt, Target, Pt_mask)
        sumcount.append(count)
        sumcover.append(cover)
    sumcount = np.array(sumcount)
    sumcount = np.sum(sumcount, axis=0)
    divide = total_step/file_step
    sumcount = sumcount/divide
    new_count = comm.gather(sumcount, root = 0)
    #print(new_count)
    #count = count.tolist()
    new_count = np.array(new_count,dtype='float')
    sumcover = np.array(sumcover)
    sumcover = np.sum(sumcover, axis=0)
    sumcover = sumcover/divide
    new_cover = comm.gather(sumcover, root = 0)
    #cover = cover.tolist()
    if rank == 0:
        print('Pt_mask Number_of_target_atom coverage(%)')
        new_count = np.array(new_count,dtype='float')
        new_count = np.sum(new_count,axis = 0)
        new_cover = np.array(new_cover)
        new_cover = np.sum(new_cover,axis = 0)
        for idx, c in enumerate(new_count):
            total_Pt = len(Pt.particles['id'][Pt.particles['mask']==idx+Pt_mask_start])
            #covered_Pt = np.sum(Pt.particles['flag'][Pt.particles['mol']==idx+1])
            print("{} {:3} {:.3}%".format(idx+Pt_mask_start, c, new_cover[idx]/total_Pt*100))

    # for debug trimmed_Pt
    #ss1 = mmps.Stream()
    #ss1.sdat = Pt
    #ss1.output_file("show_dump", 'dumppos', ["id","type","pos",'mask'] )
    #ss.output_file("trimmed_input", 'input')

if __name__ == '__main__':
    main()
