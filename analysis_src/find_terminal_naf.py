#!/usr/bin/env python3
import numpy as np
import sys
from MolCop import mmpystream as mmps
from MolCop.analysis import topology
import copy

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
        
        m_n = list(set(atoms.particles["molnum"]))
        m_n.sort()
        print("List of molecule length")
        print(m_n)

        p_num = sum((atoms.particles['molnum'] >= (target_len - _LowerCutOff)) & (atoms.particles['molnum'] <= (target_len + _UpperCutOff)))
        print("Number of Target particles is ...")
        print(p_num)
        print("Number of Target molecules is ...")
        print(p_num/target_len)

def find_terminal(atoms : mmps.Stream().sdat,start_idx=1, end_idx=362):
    atoms.create_connect_list()
    m_list = topology.create_molecule(atoms)
    atoms.add_particles_property("Nafion_endflag", _dtype=int, dim=1)
    atoms.add_particles_property("molnum", _dtype=int, dim=1)
    for l in m_list:
        num_mol = len(l)
        for ind in l:
            atoms.particles["molnum"][ind]=num_mol
    nafion_flag = (atoms.particles["molnum"] > 400) & (atoms.particles["molnum"] < 500)
    Nafion = copy.deepcopy(atoms)
    Nafion.flagconnect = True
    Nafion.trimming_particles(nafion_flag,reindex = True)
    topology.create_molecule(Nafion)
    connect_list = Nafion.connect_list
    Nafion.add_particles_property("num_of_terminalF", _dtype=int, dim=1)
    for idx,(t,c_1d) in enumerate(zip(Nafion.particles['type'],connect_list)):
        if t == 1:
            connect_count=0
            for c in c_1d:
                connect = Nafion.particles['type'][c]
                if connect == 6: #if C connect with F
                    connect_count += 1
            Nafion.particles['num_of_terminalF'][idx] = connect_count

    mol_list=list(set(Nafion.particles['mol']))
    mol_list.sort()
    for i in mol_list:
        mol_flag = Nafion.particles['mol']==i
        edge_flag=Nafion.particles['num_of_terminalF'][mol_flag] == 3
        id_of_molecule=Nafion.particles['id'][mol_flag]
        edgecarbon_id = id_of_molecule[edge_flag]
        min_edgecarbon_id = min(edgecarbon_id)
        Nafion.particles['Nafion_endflag'][min_edgecarbon_id-1]=1
        max_edgecarbon_id = max(edgecarbon_id)
        Nafion.particles['Nafion_endflag'][max_edgecarbon_id-1]=1
    #start_Nafion_idx = min(ss.sdat.particles['id'][flag])
    atoms.particles['Nafion_endflag'][nafion_flag]=Nafion.particles['Nafion_endflag']
    return atoms

if __name__ == '__main__':
    ss = mmps.Stream()
    ifn = sys.argv[1]
    ss.import_file(ifn, 'dumppos')
    ifn = sys.argv[2]
    ss.import_file(ifn, 'dumpbond')
    ss.sdat=find_terminal(ss.sdat)
    ss.output_file('show_dump', 'dumppos', ['id','type','pos','velo','force','mask','Nafion_endflag'])
