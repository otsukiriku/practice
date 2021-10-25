import copy
import math
from MolCop import mmpystream as mmps
import numpy as np
#construct object
ss1 = mmps.Stream()
ss2 = mmps.Stream()

#import input name

ss1.import_file("stable_1CB.dump", 'dumppos')
ss2.import_file("w_H2O2.dump", 'dumppos')

#construct output object
ss_out = mmps.Stream()

#define element number
#ss_out.sdat.set_elem_to_type({"O":3, "H":2})
#read element number from reaxff paramerter file
ss_out.sdat.set_elem_to_type(mmps.get_elem_to_type("para.rd"))

#define output cell size
ss_out.sdat.set_cellsize([[0, 300], [0, 300], [0, 300]])

#define pos_shift list

def g_c(atoms : mmps.Stream().sdat.particles):
    x, y, z = np.hsplit(atoms['pos'],3)
    center = np.array([x.mean(), y.mean(), z.mean()])
    return center
cent = g_c(ss1.sdat.particles)
#print(cent[0] ,cent[1], cent[2])
#print(cent)

ss_out.sdat.concate_particles(ss1.sdat.particles)
#flag = (ss2.sdat.particles['pos'][0] - cent[0])**2 +(ss2.sdat.particles['pos'][1] - cent[1])**2 +(ss2.sdat.particles['pos'][2] - cent[2])**2 < 10000
dist = (ss2.sdat.particles['pos'] - cent)**2
#同じ行のものを足す
distarr=np.sum(dist, axis=1)
flag = (distarr>10000) & (distarr<22500)
ss2.sdat.trimming_particles(flag)
ss_out.sdat.concate_particles(ss2.sdat.particles)


ss_out.output_file('newinput', 'input')
ss_out.output_file('show_dump', 'dumppos', ['id','type', 'pos'])

