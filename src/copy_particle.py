import copy
import numpy as np

from MolCop import mmpystream as mmps

#construct object
ss = mmps.Stream()

#import input name
ifn = "input"

#import input file
ss.import_file(ifn, 'input')

#construct output object
ss_out = mmps.Stream()

#define element number
#ss_out.sdat.set_elem_to_type({"O":3, "H":2})
#read element number from reaxff paramerter file
ss_out.sdat.set_elem_to_type(mmps.get_elem_to_type("para.rd"))

#define output cell size
ss_out.sdat.set_cellsize([[0, 476.413], [0, 517.263], [0,970.791]])

#define pos_shift list
pos_shift = [
    [79.700421, 76.319944, -20.2615],
    [195.959008, 80.136479, 145.684137],
    [145.93113,  223.672496, 285.1845],
    [64.8505, 153.0364,456.2845],
    [16.404836, 260.899369, 620.980014],
    [83.143586, 96.829214, 720.9855],
]
   

#ss.sdat.trimming_particles( ss.sdat.particles['type'] == 1, True)
#ss.sdat.trimming_particles( ss.sdat.particles['type'] == 1, False)

for ind, ps in enumerate(pos_shift):
    atoms = copy.deepcopy(ss.sdat.particles)
    atoms['pos'] = np.add(atoms['pos'], ps) 
    #atoms['pos'] += ps 
    atoms['mask'].fill(ind+1) 
    ss_out.sdat.concate_particles(atoms)


ss_out.output_file('newinput', 'input')
ss_out.output_file('show_dump', 'dumppos', ['id','type', 'pos'])
