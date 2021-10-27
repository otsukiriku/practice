#!/usr/bin/env python3
#trimming
from MolCop import mmpystream as mmps
from MolCop.analysis import topology

import numpy as np
import sys
import copy


ss = mmps.Stream()
ifn = sys.argv[1]
ss.import_file(ifn, 'input')
print(f"Read input file {ifn}")
ss.import_file("dump.pos.20000", 'dumppos')
print(f"Read dump file ")

ss.sdat.wrap_particles()
ss.sdat.shift_particles()

ss.import_file("dump.bond.20000", 'dumpbond')
print(f"Read bond file ")

ss.sdat.create_connect_list(0.3)

connect_list = ss.sdat.connect_list
m_list = topology.create_molecule(ss.sdat)

flag = ss.sdat.particles["mol"] == 1



ss.sdat.trimming_particles(flag, reindex=True)

O_num = np.count_nonzero(ss.sdat.particles['type']==3)
H_num = np.count_nonzero(ss.sdat.particles['type']==2)

print(f'OH:{O_num}, H:{H_num-O_num}')

ss.output_file("trimmed_dump", 'dumppos', ["id","type","pos", "mol"] )
ss.output_file("trimmed_input", 'input')

