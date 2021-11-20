import copy
import numpy as np
import sys
from MolCop import mmpystream as mmps

#construct object
ss = mmps.Stream()

#import input name


#import input file
ifn = sys.argv[1]
ss.import_file(ifn, 'dumppos')

#read element number from reaxff paramerter file
#ss.sdat.set_elem_to_type(mmps.get_elem_to_type("para.rd"))

ss.sdat.wrap_particles()
ofn = 'wrapped'+ifn
ss.output_file(ofn, 'dumppos', ['id','type', 'pos','velo'])
#ss_out.output_file('show_dump', 'dumppos', ['id','type', 'pos'])
