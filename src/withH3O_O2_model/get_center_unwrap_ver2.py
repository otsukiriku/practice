#!/usr/bin/env python3

#get_center

from MolCop import mmpystream as mmps
from MolCop.analysis import topology
import numpy as np
import sys
import copy
CB_num=1
CB_radius = 75
ss = mmps.Stream()
ss_out = mmps.Stream()
sys.path.append("/nfshome12/rotsuki/molcop/src/")
from modeling.unwrap_particles import unwrap_p

def g_c(atoms : mmps.Stream().sdat, flag=None):
    if flag is None:
        center = np.array(atoms.particles['pos'].mean(axis=0))
    else:
        center = np.array(atoms.particles['pos'][flag].mean(axis=0))
    return center

def correct_center(atoms : mmps.Stream().sdat, center:np.array):
    half_cell=CB_radius +80
    for idx,cent in enumerate(center):
        flag=ss.sdat.particles['mask']==idx+1
        over_flag = (atoms.particles['pos'][flag]-cent) > half_cell
        under_flag = (atoms.particles['pos'][flag]-cent) < -half_cell
        atoms.particles['pos'][flag] -= over_flag * atoms.newcell
        atoms.particles['pos'][flag] += under_flag * atoms.newcell

if __name__ == '__main__':
    print(f'{CB_num} Partcles Radius:  110.000000 A')
    args = sys.argv
    ifn = args[1]
    ss.import_file(ifn, 'dumppos')
    #ss.sdat.create_connect_list()
    tmpflag = ss.sdat.particles["id"]< 122296 
    tmpflag &= ss.sdat.particles["type"]==1
    ss.sdat.particles["mask"][tmpflag] = 1
    center_list=[]
    for mask in range(1,CB_num+1):
        fl=ss.sdat.particles['mask']==mask
        c_out = g_c(ss.sdat,fl)
        center_list.append(c_out)
        #print(mask-1, end = " ")
        #print(*c_out)
    center_list=np.array(center_list)
    correct_center(ss.sdat,center_list)

    for mask in range(1,CB_num+1):
        fl=ss.sdat.particles['mask']==mask
        c_out = g_c(ss.sdat,fl)
        print(mask-1, end = " ")
        print(*c_out)

    ss.output_file("show_dump","dumppos", ["id","type", "pos","mask"])
