#!/usr/bin/env python3
#get_center
from MolCop import mmpystream as mmps
import numpy as np
import sys
import copy
ss = mmps.Stream()
def g_c(atoms : mmps.Stream().sdat.particles):
        x, y, z = np.hsplit(atoms['pos'],3)
        center = np.array([x.mean(), y.mean(), z.mean()])
        return center
if __name__ == '__main__':
    args = sys.argv
    ifn = args[1]
    ss.import_file(ifn, 'input')
    dcell=ss.sdat.cell[2]
    pos_shift = [0, 0, -dcell]
    print(pos_shift)
    for mask in range(1,7):
        sdat_mask = copy.deepcopy(ss.sdat)
        sdat_mask.trimming_particles(sdat_mask.particles['mask']==mask, False)

        flag = ((sdat_mask.particles['mask']==1) & (sdat_mask.particles['pos'][:,2]>750))
        #print(flag)
        sdat_mask.particles['pos'][flag]+=pos_shift

        #sdat_mask['pos'][fl] = np.add(sdat_mask['pos'][fl],pos_shift)
        #sdat_mask.
        #sdat_mask.shift_particles(shift=pos_shift)
        atoms = sdat_mask.particles
        c_out = g_c(atoms)
        print(mask - 1, end=' ')
        [print(co, end=' ') for co in c_out]
        print(end='\n')
        del sdat_mask
