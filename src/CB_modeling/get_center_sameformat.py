#!/usr/bin/env python3

#get_center

from MolCop import mmpystream as mmps
import numpy as np
import sys
import copy
CB_num=4
ss = mmps.Stream()

def g_c(atoms : mmps.Stream().sdat, flag=None):
    if flag is None:
        center = np.array(atoms.particles['pos'].mean(axis=0))
    else:
        center = np.array(atoms.particles['pos'][flag].mean(axis=0))
    return center

if __name__ == '__main__':
    args = sys.argv
    ifn = args[1]
    ss.import_file(ifn, 'input')
    obj=ss.sdat
    cell_len=obj.cell[2]
    pos_shift = [0,0,cell_len]
    CB = ss.sdat.particles
    flag = ((CB['mask']==1) & (CB['pos'][:,2]>cell_len*2/3))
    CB['pos'][flag]-=pos_shift
    flag = ((CB['mask']==CB_num) & (CB['pos'][:,2]<cell_len*1/3))
    CB['pos'][flag]+=pos_shift
    print(f'{CB_num} Partcles Radius:  110.000000 A')
    for mask in range(1,CB_num+1):
        flag = ss.sdat.particles['mask'] == mask
        c_out = g_c(ss.sdat, flag)
        print(mask-1, end = " ")
        print(*c_out)
    if len(args) >= 3:
        ofn = args[2]
        li = c_out.tolist()
        with open(ofn, mode='w') as f:
            for q in li:
                f.write(str(q)+' ')
            f.write('\n')
