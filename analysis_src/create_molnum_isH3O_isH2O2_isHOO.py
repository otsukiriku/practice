#!/usr/bin/env python3
from MolCop import mmpystream as mmps
from MolCop.analysis import topology
from E_T.func import molecule
from E_T.func import read_center as rc
from E_T.func import get_center as gc
import numpy as np
import sys
import copy
import math
sys.path.append("/nfshome12/rotsuki/molcop/src/")
#print(sys.path)
import modeling.unwrap_particles as u_p
import evaluate_structure.get_center as g_c

#for MPI
import os

#********************************
Pt_cluster_num=309
Pt_num=14
Pt_mask_start = 7
#cut off distance from Pt surface
cutoff = 9
target_type = 6
# 5:SofNafion 6:FofNadion
#********************************

#******************************
target_name = "H3O"
search_target = ["H","O","H","H"]
target_name2 = "H2O2"
search_target2 = ["H","O","O","H"]
target_name3 = "HOO"
search_target3 = ["H","O","O"]
#******************************

def find_molecule(atoms:mmps.Stream().sdat,target_name,search_target):
    search_target = np.array(sorted(search_target))
    atoms.type_to_elem.clear()
    atoms.type_to_elem = {1:"C",2:"H",3:"O",4:"Pt",5:"S",6:"F",12:"H",13:"O",14:"H",15:"O",16:"O",21:"C"}

    if len(atoms.connect_list) == 0:
        atoms.create_connect_list()
    if "molnum" not in atoms.particles:
        molecule.create_molnum(atoms)
    if "mol" not in atoms.particles:
        topology.create_molecule(atoms)

    target_name = "is" + target_name
    if target_name not in atoms.particles:
        atoms.add_particles_property(target_name,int,dim=1)
    molnum_flag = atoms.particles["molnum"] == len(search_target)
    target_type = atoms.particles["type"][molnum_flag]
    target_mol = atoms.particles["mol"][molnum_flag]
    new_property = atoms.particles[target_name][molnum_flag]
    target_mol_set = set(target_mol)

    for t_m in target_mol_set:
        flag = target_mol == t_m
        tmp_target_type = target_type[flag]
        judge_lst = [ atoms.type_to_elem[t] for t in tmp_target_type ]
        judge_lst = sorted(judge_lst)
        judge_lst = np.array(sorted(judge_lst))
        if (search_target == judge_lst).all():
            new_property[flag] = 1
    atoms.particles[target_name][molnum_flag] = new_property

    
ss = mmps.Stream()
ss.import_file('../config.rd', 'config')

start_step = int(sys.argv[1])
end_step = int(sys.argv[2])
total_step = end_step - start_step

timestep = []
isH3O = []
isH2O2 = []
isHOO = []

for step in range(start_step, end_step+1, ss.sdat.file_step):
#for step in range(start_step, end_step+1, 4000):
    ss = mmps.Stream()

    #timestep.append(step)
    #isH3O.append(np.sum(ss.sdat.particles["isH3O"]))
    #isH2O2.append(np.sum(ss.sdat.particles["isH2O2"]))

    ifn = 'colored.pos.' + str(step)
    if not os.path.isfile(ifn):
        print(f"{ifn} is not exit")
        continue
    ss.import_file(ifn, 'dumppos')
    ifn1 = '../dump.bond.' + str(step)
    if not os.path.isfile(ifn1):
        print(f"{ifn1} is not exit")
        continue
    ss.import_file(ifn1, 'dumpbond')
    ss.sdat.create_connect_list()
    molecule.create_molnum(ss.sdat)
    mask = []
    for i in range(1,Pt_num+1):
        mask.append([ i for _ in range(Pt_cluster_num) ])
    np.array(mask)
    mask = np.ravel(mask)
    ss.sdat.particles["mask"][ss.sdat.particles["type"]==4] = mask

    find_molecule(ss.sdat,target_name,search_target)
    find_molecule(ss.sdat,target_name2,search_target2)
    find_molecule(ss.sdat,target_name3,search_target3)
    ofn = "withmolnum_" + ifn
    ss.output_file(ofn, "dumppos",["id","type","pos","velo","mask","molnum","is"+target_name, "is"+target_name2,"is"+target_name3])
