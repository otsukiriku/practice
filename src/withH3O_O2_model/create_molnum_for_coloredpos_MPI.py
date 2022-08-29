#!/usr/bin/env python3
from MolCop import mmpystream as mmps
from MolCop.analysis import topology
from E_T.func import molecule
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



def main():
    # for MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    if rank == 0:
        print(f"Input process : {size}")
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
        ifn = 'colored.pos.' + str(step)
        if not os.path.isfile(ifn):
            print(f"{ifn} is not exit")
            continue
        ss.import_file(ifn, 'dumppos')
        ifn1 = '../dump.bond.' + str(step)
        if not os.path.isfile(ifn1):
            print(f"{ifn1} is not exit")
            continue
        ss.import_file(ifn1, 'dumpbond')
        ss.sdat.create_connect_list()
        molecule.create_molnum(ss.sdat)
        ofn = "withmolnum_" + ifn
        ss.output_file(ofn, "dumppos",["id","type","pos","velo","mask","molnum"])

if __name__ == '__main__':
    main()
