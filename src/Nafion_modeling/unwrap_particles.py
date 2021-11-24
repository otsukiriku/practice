
import copy
import numpy as np
import sys
import os.path
from MolCop import mmpystream as mmps
from MolCop.analysis import topology 

#sys.path.append("/nfshome12/rotsuki/molcop/src/evaluate_structure")
sys.path.append("~/molcop/src/evaluate_structure")


print(sys.path)
#from get_center import g_c
import get_center as g_c


#construct object
ss = mmps.Stream()

#import input name


#import input file
ss.import_file("dump.pos.0", 'dumppos')
ss.import_file("dump.bond.0", 'dumpbond')
ss.sdat.create_connect_list(0.3)
topology.unwrap_molecule(ss.sdat)
m_list = topology.create_molecule(ss.sdat)

atoms = ss.sdat.particles

for ind, m in enumerate(m_list):
    flag = atoms["mol"] == ind+1
    gc = g_c.g_c(atoms, flag)
    print(gc)

ss.output_file("unwrap_dump.pos.0", "dumppos", ["id", "type", "pos"])


#read element number from reaxff paramerter file
#ss.sdat.set_elem_to_type(mmps.get_elem_to_type("para.rd"))

