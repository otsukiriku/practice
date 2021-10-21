#!/usr/bin/env python3
#trimming
from MolCop import mmpystream as mmps
from MolCop.analysis import topology

import numpy as np
import sys
import copy



ss = mmps.Stream()
#ifn = sys.args[1]
ss.import_file("dump.pos.0", 'dumppos')
ss.import_file("dump.bond.0", 'dumpbond')

ss.sdat.create_connect_list(0.3)

connect_list = ss.sdat.connect_list
m_list = topology.create_molecule(ss.sdat)

flag = ss.sdat.particles["mol"] == 1



ss.sdat.trimming_particles(flag, reindex=True)

ss.output_file("trimmed_dump", 'dumppos', ["id","type","pos", "mol"] )
ss.output_file("trimmed_input", 'input')

