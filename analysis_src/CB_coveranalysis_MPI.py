#!/usr/bin/env python3
import os
import sys
import copy
from mpi4py import MPI
from MolCop import mmpystream as mmps
from MolCop.analysis import topology, neighbor
import numpy as np
import time

#****************************************
Pt_cluster_num=309
Pt_num=56
CB_maskstart = 7
#cut off distance from Pt surface
cutoff = 6
target_type = 6
# 5:SofNafion 6:FofNadion
#****************************************


def calc_coverage(target_pos:np.ndarray ,atoms:mmps.Stream().sdat):
    count=0
    cover=0
    cover_flag = []
    limit_calcarea
    for p in atoms.particles['pos']:
        dist = (target_pos - p)**2
        dist = np.sum(dist, axis=1)
        #[print(math.sqrt(d)) for d in dist if d < cutoff**2]
        flag = dist < (cutoff**2)
        sumflag = np.sum(flag)
        coverflag=0
        cfl = 0
        if sumflag != 0:
            coverflag = 1
            cfl = 1
        cover_flag.append(cfl)
        cover+=coverflag
        count+=sumflag
    print(cover_flag)
    return count,cover,cover_flag

def main():

    #Init MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    if rank == 0:
        print(f"Input process : {size}")

    ss = mmps.Stream()
    ss.import_file('config.rd', 'config')

    start_step = int(sys.argv[1])
    end_step = int(sys.argv[2])
    total_step = end_step - start_step

    file_step = ss.sdat.file_step
    dstep = int(total_step/(file_step*size)) + 1

    init_step = dstep * rank * file_step + start_step
    final_step = min(dstep * (rank+1) * file_step + start_step, total_step + start_step)

    if final_step == total_step:
        final_step += 1

    for step in range(init_step, final_step, ss.sdat.file_step):

        ifn_pos = '../colored.pos.' +  str(step)
        if not os.path.isfile(ifn_pos):
            print(f'{ifn_pos} is not exit')
            continue
        ss.import_file(ifn_pos, 'dumppos')

        """
        #for CB
        CB = copy.deepcopy(ss.sdat)
        startptidx = np.min(ss.sdat.particles['id'][ss.sdat.particles['type']==4])
        CB_flag = ((ss.sdat.particles['type']==21)|(ss.sdat.particles['type']==22))  & (ss.sdat.particles['id']<startptidx)
        CB.flagconnect = True
        CB.trimming_particles(CB_flag, reindex=True)
        # for target
        target_flag = ss.sdat.particles['type'] == target_type
        tar_pos = copy.deepcopy(ss.sdat.particles['pos'][target_flag])
        ss.sdat.add_particles_property('coverflag', _dtype=int, dim=1)
        ss.sdat.particles['coverflag'][CB_flag] = cover_flag
        ofs = "CBcover_colored.pos."+str(step)
        ss.output_file(ofs, 'dumppos',["id","type","pos","velo","force","q","mask", "coverflag"])
        """
        startptidx = np.min(ss.sdat.particles['id'][ss.sdat.particles['type']==4])
        CB_flag = ((ss.sdat.particles['type']==21)|(ss.sdat.particles['type']==22))  & (ss.sdat.particles['id']<startptidx)
        target_flag = ss.sdat.particles['type'] == target_type
        neighbor_flag = CB_flag | target_flag
        target = copy.deepcopy(ss.sdat)
        target.trimming_particles(neighbor_flag,reindex=True)
        n_list = neighbor.make_neighbor(target, 6)
        cover_flag = []
        for idx,n_1d in enumerate(n_list):
            if idx == startptidx-1:
                break
            for n in n_1d:
                if n >= startptidx - 1:
                    flag = True
                    break
                else:
                    flag = False
            if flag == False:
                cover_flag.append(0)
            else:
                cover_flag.append(1)

        #print("CB_particle_num:{}".format(startptidx))
        #print("length of cover_flag:{}".format(len(cover_flag)))
        cover_flag=np.array(cover_flag)
        if 'coverflag' not in ss.sdat.particles:
            ss.sdat.add_particles_property('coverflag', _dtype=int, dim=1)
        ss.sdat.particles['coverflag'][CB_flag] = cover_flag
        ofs = "CBcover_colored.pos."+str(step)
        ss.output_file(ofs, 'dumppos',["id","type","pos","velo","force","q","mask", "coverflag"])
        
if __name__ == '__main__':
    a=time.time()
    main()
    b=time.time()
    print("time{}".format((b-a)))

