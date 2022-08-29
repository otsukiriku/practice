
import sys
import numpy as np
from MolCop import mmpystream as mmps


#construct object
ss = mmps.Stream()

if len(sys.argv) <= 1:
    print("No arguments are detected")
    sys.exit()


#import dumppos file
ifn = sys.argv[1]

ss.import_file(ifn, 'input')

new_id=[i+1 for i in range(ss.sdat.total_particle)]
ss.sdat.particles['id'] = np.array(new_id)
ss.sdat.sort_particles_by_id()

ofn="new_"+ifn
ss.output_file(ofn,'input')
#ss.output_file(ofn,'dumppos',["id","type","pos","velo","mask"])
