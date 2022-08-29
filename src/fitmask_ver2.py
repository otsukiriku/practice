import sys
import copy
import numpy as np
from MolCop import mmpystream as mmps

#construct object
ss = mmps.Stream()
ss1 = mmps.Stream()
ss_out = mmps.Stream()

if len(sys.argv) <= 1:
    print("No arguments are detected")
    sys.exit()

#import input files
if len(sys.argv) >= 3:
    ifn = sys.argv[2]
    ss1.import_file(ifn, 'input')
    print(f"Read input file {ifn}")

#import dumppos file
ifn = sys.argv[1]
#import input file

ss.import_file(ifn, 'dumppos')
print(f"Read file {ifn}")


ss_len=len(ss.sdat.particles["id"])
ss1_len=len(ss1.sdat.particles["id"])


if ss_len>ss1_len:
    differ=ss_len-ss1_len
    z=np.zeros(differ, dtype=int)
    inp_mask=ss1.sdat.particles["mask"]
    inp_mask=np.append(inp_mask, z)
    ss.sdat.particles["mask"]=inp_mask
    ss_out.sdat = copy.deepcopy(ss.sdat)

if ss_len<ss1_len:
    differ=ss1_len-ss_len
    z=np.zeros(differ, dtype=int)
    inp_mask=ss.sdat.particles["mask"]
    inp_mask=np.append(inp_mask, z)
    ss1.sdat.particles["mask"]=inp_mask
    ss_out.sdat = copy.deepcopy(ss1.sdat)

ss_out.sdat.set_elem_to_type(mmps.get_elem_to_type("para.rd"))

ofs = "newinput_from_" + ifn
print(f"Read file {ifn}")

ofs = "newinput_from_" + ifn
ss_out.output_file(ofs, 'input')
