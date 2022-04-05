import copy
import math
import sys
from MolCop import mmpystream as mmps
from MolCop.analysis import topology
from E_T.func import molecule
from E_T.io import read_center as rc
import numpy as np
#construct object
ss1 = mmps.Stream()
ss2 = mmps.Stream()

#********************************
CB_num=10
Pt_num=57
water_maxradi=83.5
water_minradi=76
CB_radi=75
#********************************

#import input name

ss1.import_file(sys.argv[1], 'dumppos')
ss1.import_file(sys.argv[2], 'dumpbond')
ss2.import_file("inside_water.dump", 'dumppos')
#ss2.import_file("inside_water.bond", 'dumpbond')

#construct output object
ss_out = mmps.Stream()

#read element number from reaxff paramerter file
ss_out.sdat.set_elem_to_type(mmps.get_elem_to_type("para.rd"))

#define output cell size
cell_len=ss1.sdat.cell
pos_shift = [0,0,cell_len[2]]
ss_out.sdat.cell=ss1.sdat.cell
ss_out.sdat.dcell=ss1.sdat.dcell
ss_out.sdat.newcell=ss1.sdat.newcell
ss_out.sdat.concate_particles(ss1.sdat.particles)

def g_c(atoms : mmps.Stream().sdat.particles, flag=None):
    if flag is None:
        center = np.array([atoms['pos'].mean(axis=0)])
    else:
        center = np.array([atoms['pos'][flag].mean(axis=0)])
    return center

def g_c_by_mask(atoms : mmps.Stream().sdat.particles, mask):
    flag = atoms["mask"] == mask
    center = atoms["pos"][flag].mean(axis=0)
    return center

CB_G = rc.r_c()

#extend CB_G
CB_G = np.array(CB_G)
CB_org_G = np.array(CB_G)
cell= np.array(ss1.sdat.cell)
#print("CB_G")
#print(CB_G)
wrap_list=np.array([[1,0,0],[0,1,0],[0,0,1]])
for lst in wrap_list:
    p_shift=CB_org_G+cell*lst
    m_shift=CB_org_G-cell*lst
    CB_G=np.append(CB_G,p_shift,axis=0)
    CB_G=np.append(CB_G,m_shift,axis=0)
    #if ( (cell+CB_radi + 20) > p_shift).all():
    #    CB_G=np.append(CB_G,p_shift,axis=0)
    #if ( (-CB_radi -20) < m_shift).all():
    #    CB_G=np.append(CB_G,m_shift,axis=0)
print("CB_G")
print(CB_G)

ss2.sdat.shift_particles()
H2O = ss2.sdat
#delete particles from center
ss2.sdat.replicate_particles([3,3,3])
#ss2.sdat.cell=ss1.sdat.cell
#ss2.sdat.dcell=ss1.sdat.dcell
#ss2.sdat.newcell=ss1.sdat.newcell

over_flag = (ss2.sdat.particles['pos'] > cell)
over_flag=np.sum(over_flag,axis=1)
over_flag.flatten()
over_flag = over_flag == 0
ss2.sdat.trimming_particles(over_flag, reindex=True)

elem_num=len(ss2.sdat.particles['id'])
flag_2d = np.zeros((1,elem_num), dtype=int)
for c in CB_G:
    dist = (H2O.particles['pos'] - c)**2
    #同じ行のものを足す
    distarr=np.sum(dist, axis=1)
    flag = ((water_minradi**2<distarr) & (distarr<water_maxradi**2))
    #bool -> int 
    flag = flag.astype(int)
    flag_insideCB = (distarr<(CB_radi**2+500))
    print(flag_insideCB)
    flag =  flag + (-1)*flag_insideCB.astype(int)  
    flag_2d = np.append(flag_2d, [flag], axis=0)
#print(flag_2d) 
flag = np.sum(flag_2d, axis=0) > 0
#print(flag)
ss2.sdat.trimming_particles(flag)

"""
#trim_inside water_after wrapping
elem_num=len(ss2.sdat.particles['id'])
flag_2d = np.zeros((1,elem_num), dtype=bool)
#print(flag)
for c in CB_G:
    dist = (ss2.sdat.particles['pos'] - c)**2
    distarr=np.sum(dist, axis=1)
    flag_insideCB = (distarr<(CB_radi**2+500))
    flag_2d|=flag_insideCB
    #flag_2d = np.append(flag_2d, [flag], axis=0)
print(flag_2d)
flag = np.logical_not(flag_2d[0])
#print(flag)
#print(ss2.sdat.particles['id'].shape)
ss2.sdat.trimming_particles(flag)
"""

#trim over cell
#flag = (0<H2O.particles['pos'][:,2]) & (H2O.particles['pos'][:,2]<cell_len[2])   
#ss2.sdat.trimming_particles(flag)

#trim over cell
#Pt_flag = ss1.sdat.particles['type']==4
#Pt_pos = ss1.sdat.particles['pos'][Pt_flag]

Pt_mask_start=min(ss1.sdat.particles['mask'][ss1.sdat.particles['type']==4])
Pt_list=[i for i in range(Pt_mask_start,Pt_mask_start+Pt_num)]

#this is for wrapped Pt
#flag=((ss1.sdat.particles['mask']==10) & (ss1.sdat.particles['pos'][:,2]>cell_len[2]*2/3))
#ss1.sdat.particles['pos'][flag]-=pos_shift

#ss1.sdat.create_connect_list()
#topology.unwrap_molecule(ss1.sdat)
Pt_G = np.array([g_c_by_mask(ss1.sdat.particles, l) for l in Pt_list])
print(Pt_G)

Pt_G = rc.r_c(filename="Pt_center.txt")
print(Pt_G)

flag_2d_Pt = np.zeros((1,len(ss2.sdat.particles['id'])), dtype=int)
for pt_g in Pt_G:
    dist = (ss2.sdat.particles['pos']-pt_g)**2
    distarr=np.sum(dist, axis=1)
    dist_flag2 = distarr > 14*14
    flag_2d_Pt = np.append(flag_2d_Pt, [dist_flag2], axis=0)
dist_flag2 = np.sum(flag_2d_Pt.astype(int), axis=0)
max_value=dist_flag2.max()
dist_flag2 = (dist_flag2 == max_value)
ss2.sdat.trimming_particles(dist_flag2, reindex=True)


"""
flag_2d_Pt = np.zeros((1,len(ss2.sdat.particles['id'])), dtype=int)
for pt_g in Pt_G:
    dist = (ss2.sdat.particles['pos']-pt_g)**2
    distarr=np.sum(dist, axis=1)
    dist_flag2 = distarr > 14*14
    flag_2d_Pt = np.append(flag_2d_Pt, [dist_flag2], axis=0)
dist_flag2 = np.sum(flag_2d_Pt.astype(int), axis=0)
max_value=dist_flag2.max()
dist_flag2 = (dist_flag2 == max_value)
ss2.sdat.trimming_particles(dist_flag2, reindex=True)
"""

#trim close Pt 
ss2.sdat.add_particles_property('mask', int, 1)

#ss2.sdat.create_connect_list()
ss_out.sdat.concate_particles(ss2.sdat.particles)
ss_out.sdat.wrap_particles()
print("number of water particles")
print(len(ss2.sdat.particles['id'])/3)
ss_out.output_file('test_inside_water.rd', 'input')
ss_out.output_file('show_dump', 'dumppos', ['id','type','pos','velo','mask'])

"""
#debug_show_CB_G
debug = mmps.Stream()
for idx, r in enumerate(CB_G):
    if "id" not in debug.sdat.particles:
        debug.sdat.add_particles_property("id", _dtype=int)
    debug.sdat.particles['id']=np.append(debug.sdat.particles['id'],idx)
    if "pos" not in debug.sdat.particles:
        debug.sdat.add_particles_property("pos",dim=3)
    debug.sdat.particles['pos']=np.append(debug.sdat.particles['pos'],[r],axis=0)
    if "type" not in debug.sdat.particles:
        debug.sdat.add_particles_property("type", _dtype=int)
    debug.sdat.particles['type']=np.append(debug.sdat.particles['type'],99)

ss_out.sdat.concate_particles(debug.sdat.particles)
ss_out.output_file('debug_dump', 'dumppos', ['id','type','pos'])
"""
