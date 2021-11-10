import sys

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
else if dump_len>input_len:
    differ=dump_len-input_len
    inp_mask=ss1.sdat.particles["mask"]
    inp_mask=[np.append(inp_mask, 0) for i in differ]
else:
    print("input has larger particle number!!!")i
    quit()


ofs = "newinput_from_" + ifn
print(f"Read file {ifn}")

ofs = "newinput_from_" + ifn
ss.output_file(ofs, 'input')
