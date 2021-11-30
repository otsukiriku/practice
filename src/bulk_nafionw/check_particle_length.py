#!/usr/bin/env python3
import numpy as np
import sys
from MolCop import mmpystream as mmps
from MolCop.analysis import topology

def p_len(atoms : mmps.Stream().sdat, target_len=None, _UpperCutOff = 0 , _LowerCutOff = 0):
    
    if target_len is None:
        print("You need to set Target length of particle!")
        exit()
    else:
        atoms.create_connect_list(0.3)
        connect_list = atoms.connect_list
        m_list = topology.create_molecule(atoms)

        atoms.add_particles_property("molnum", _dtype=int, dim=1)
        for l in m_list:
            num_mol = len(l)
            for ind in l:
                atoms.particles["molnum"][ind]=num_mol
        
        print("List of molecule length")
        m_n = list(set(atoms.particles["molnum"]))
        m_n.sort()
        print(m_n)

        mn_list = []
        if _UpperCutOff != 0 | _LowerCutOff != 0:
            for i in range(target_len - _LowerCutOff , target_len + _UpperCutOff +1):
                mn_list.append(i)
        else:
            mn_list = [target_len]

        flag_2d = np.zeros((1,len(atoms.particles['id'])), dtype=int)
        for i in mn_list:
            flag = atoms.particles["molnum"] == i
            #print(flag)
            flag = flag.astype(int)
            flag_2d = np.append(flag_2d, [flag], axis=0)
        flag1d = np.sum(flag_2d, axis=0) != 0
        #print(flag1d)

        print("Number of Target particles is ...")
        print(len(atoms.particles["id"][flag1d])/target_len)
        
if __name__ == '__main__':
    ss = mmps.Stream()
    ifn = sys.argv[1]
    ss.import_file(ifn, 'dumppos')
    ifn = sys.argv[2]
    ss.import_file(ifn, 'dumpbond')
    #example target molecule is Nafion (Molecule is composed of 416 particles)
    t_len = 416
    #Nafion has 6 sidechains, Assuming all sidechain coordinate with water, 3×6 = 18
    upper = 18
    #Nafion has 6 sidechains, Assuming H particle away from all side chains 1×6 = 6
    lower = 6
    p_len(ss.sdat ,target_len=t_len, _UpperCutOff=upper, _LowerCutOff=lower)
    
