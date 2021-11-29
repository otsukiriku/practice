#!/usr/bin/env python3
import numpy as np
import copy
from MolCop import mmpystream as mmps
from MolCop.analysis import topology
end_of_nafion = 89856

ss = mmps.Stream()
ss_out = mmps.Stream()

ss.import_file('dump.pos.0', 'dumppos')
ss.import_file('dump.bond.0', 'dumpbond')

ss.sdat.set_elem_to_type(mmps.get_elem_to_type("./para.rd"))

ss.sdat.sort_particles_by_id()
ss_out.sdat.cell = ss.sdat.cell
ss_out.sdat.dcell = ss.sdat.dcell
ss_out.sdat.newcell = ss.sdat.newcell
flag=ss.sdat.particles['id'] > end_of_water
print(ss.sdat.particles['id'][flag])
water = copy.deepcopy(ss.sdat)
water.trimming_particles(flag)
ss_out.sdat.concate_particles(water.particles)

flag=ss.sdat.particles['id'] < (end_of_water+1)
ss.sdat.trimming_particles(flag)
ss_out.sdat.concate_particles(ss.sdat.particles)

ss_out.output_file('newnafionw.dump', 'dumppos', ['id','type','pos','velo'])
#ss_out.output_file("withwaterlayer.rd", 'input')


