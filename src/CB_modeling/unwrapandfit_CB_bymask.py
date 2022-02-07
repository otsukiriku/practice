
import copy
import numpy as np
import sys
import os

from MolCop import mmpystream as mmps
from MolCop.analysis import topology
#sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append("/nfshome12/rotsuki/molcop/src/")
#print(sys.path)

#from get_center import g_c
import evaluate_structure.get_center as g_c

def unwrap_p(atoms : mmps.Stream().sdat):
    atoms.create_connect_list(0.3)
    topology.unwrap_molecule(atoms)
    m_list = topology.create_molecule(atoms)

    center_list=[]
    for ind, m in enumerate(m_list):
        flag = atoms.particles["mol"] == ind+1
        gc = g_c.g_c(atoms.particles, flag)
        #np.appendは遅いので，いったんpythonのリストにしてから再変換する
        gc = gc.tolist()
        center_list.append(gc[0])
    center_list = np.asarray(center_list)
    #print(center_list)

    for ind,c_list in enumerate(center_list):
        flag = atoms.particles["mol"] == ind+1
        over_cell = c_list >= atoms.cell
        under_cell = c_list < atoms.dcell
        atoms.particles['pos'][flag] -= over_cell * atoms.newcell
        atoms.particles['pos'][flag] += under_cell * atoms.newcell
    #print(over_list)

if __name__ == '__main__':
    ss = mmps.Stream()
    ss_wrapped = mmps.Stream()
    ifn = sys.argv[1]
    ss.import_file(ifn, 'dumppos')
    ifn = sys.argv[2]
    ss.import_file(ifn, 'dumpbond')
    mask_list=list(set(ss.sdat.particles['mask']))
    
    #unwrap_CB
    cell=np.array(ss.sdat.cell)
    for i in mask_list:
        flag=(ss.sdat.particles['mask']==i)
        center=g_c.g_c(ss.sdat.particles, flag)
        center = center[0]
        DistFromCent = ss.sdat.particles['pos'][flag] - center
        shift_p = DistFromCent > cell*0.5
        shift_m = DistFromCent < -cell*0.5
        ss.sdat.particles['pos'][flag] -= shift_p*cell
        ss.sdat.particles['pos'][flag] += shift_m*cell
    #center.append(*g_c.g_c(ss.sdat.particles, flag))
    #center=np.array(center)

    #for i in mask_list:

    max_arr=[]
    min_arr=[]
    for i in range(0,3):
        arr = ss.sdat.particles['pos'][:,i:i+1]
        arr=np.ravel(arr)
        max_arr.append(np.max(arr, axis=0)+10)
        min_arr.append(np.min(arr, axis=0)-10)
    ss.sdat.cell = max_arr
    ss.sdat.dcell = min_arr
    #atoms.set_elem_to_type(mmps.get_elem_to_type("para.rd"))
    ss.sdat.shift_particles()
    ss.output_file("unwrapped_pos", 'dumppos', ["id","type","pos","velo"] )

