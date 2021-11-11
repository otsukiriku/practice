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
    ss.import_file(ifn, 'newinput')
    print("x, y, z")
    for mask in range(1,7):
        sdat_mask = copy.deepcopy(ss.sdat)
        sdat_mask.trimming_particles(sdat_mask.particles['mask']==mask, False)
        atoms = sdat_mask.particles
        c_out = g_c(atoms)
        print(c_out)
        del sdat_mask
