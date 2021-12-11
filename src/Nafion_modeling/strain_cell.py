
import copy
import numpy as np
import sys
import os

from MolCop import mmpystream as mmps
from MolCop.analysis import topology
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append("/nfshome12/rotsuki/molcop/src/")
#print(sys.path)
    
ss = mmps.Stream()
ifn = sys.argv[1]
ss.import_file(ifn, 'dumppos')
#ifn = sys.argv[2]
#ss.import_file(ifn, 'dumpbond')

ss.sdat.strain_cell([0.0,0.0,0.0087251])

#atoms.set_elem_to_type(mmps.get_elem_to_type("para.rd"))
ss.output_file("strained_pos", 'dumppos', ["id","type","pos","velo"] )
ss.output_file("strained_input", 'input')
