import copy
import sys
from MolCop import mmpystream as mmps

#construct object
ss = mmps.Stream()

#import input name
args = sys.argv
ifn = args[1]
#import input file
ss.import_file(ifn, 'input')

flag = ss.sdat.particles['pos'][:,2] < ss.sdat.cell[2]
ss.sdat.trimming_particles(flag)
flag = ss.sdat.particles['pos'][:,2] > ss.sdat.dcell[2]
ss.sdat.trimming_particles(flag, reindex=True)

ss.output_file('newinput', 'input')



#construct output object
#ss_out = mmps.Stream()

#define element number
#read element number from reaxff paramerter file
#ss.sdat.set_elem_to_type(mmps.get_elem_to_type("para.rd"))

#ss.sdat.particles['pos']

ss.output_file('show_dump', 'dumppos', ['id','type', 'pos'])
