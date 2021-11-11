import sys
import numpy as np
from MolCop import mmpystream as mmps

#construct object
ss = mmps.Stream()
ss1 = mmps.Stream()

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


dump_len=len(ss.sdat.particles["mask"])
input_len=len(ss1.sdat.particles["mask"])

if dump_len == input_len:
    ss.sdat.particles["mask"]=ss1.sdat.particles["mask"]

elif dump_len>input_len:
    differ=dump_len-input_len
    z=np.zeros(differ, dtype=int)
    inp_mask=ss1.sdat.particles["mask"]
    inp_mask=np.append(inp_mask, z)
    ss.sdat.particles["mask"]=inp_mask

else:
    print("input should has larger particle number!!!")
    quit()

ss.sdat.set_elem_to_type(mmps.get_elem_to_type("para.rd"))

ofs = "newinput_from_" + ifn
print(f"Read file {ifn}")

ofs = "newinput_from_" + ifn
ss.output_file(ofs, 'input')
