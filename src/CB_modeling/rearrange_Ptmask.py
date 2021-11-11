
import sys
import numpy as np
from MolCop import mmpystream as mmps

Pt_num = 56
Pt_cluster_pnum=309

#construct object
ss = mmps.Stream()
ss1 = mmps.Stream()

if len(sys.argv) <= 1:
    print("No arguments are detected")
    sys.exit()


#import dumppos file
ifn = sys.argv[1]

ss.import_file(ifn, 'dumppos')
print(f"Read file {ifn}")

Pt=ss.sdat
Pt_list=[i for i in range(7,7+Pt_num)]
Ptmask=Pt.particles['mask'][Pt.particles['type'] == 4]
Ptmask=np.reshape(Ptmask, (Pt_num, Pt_cluster_pnum))

for idx, msk in enumerate(Ptmask):
    #一列ずつPtmaskから切り取る．
    msk[:] = Pt_list[idx]
#print(Ptmask)

Ptmask = Ptmask.flatten()
ss.sdat.particles['mask'][Pt.particles['type'] == 4] = Ptmask

ofn="new_"+ifn
ss.output_file(ofn,'dumppos',["id","type","pos","velo","mask"])
