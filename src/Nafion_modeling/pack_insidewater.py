import copy
import math
from MolCop import mmpystream as mmps
from MolCop.analysis import topology
import numpy as np
#construct object
CB_num=4
Pt_num=56
ss1 = mmps.Stream()
ss2 = mmps.Stream()

#import input name

ss1.import_file("trimmed_dump", 'dumppos')
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
ss_out.sdat.newcell=ss1.sdat.cell
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

CB_list=[i+1 for i in range(CB_num)]
CB_G = np.empty((0,3))
for mask in CB_list:
    #for wrap particles
    CB = ss1.sdat.particles
    flag = ((CB['mask']==1) & (CB['pos'][:,2]>cell_len[2]*2/3))
    CB['pos'][flag]-=pos_shift
    flag = ((CB['mask']==CB_num) & (CB['pos'][:,2]<cell_len[2]*1/3))
    CB['pos'][flag]+=pos_shift

    flag = CB["mask"] == mask
    #print(CB["mask"]) 
    #print(flag) 
    CB_g = [CB["pos"][flag].mean(axis=0)]
    CB_G = np.append(CB_G, CB_g, axis=0)
a=CB_num-1
bottom = CB_G[CB_num-1,] - np.array(pos_shift)
CB_G = np.append(CB_G, [bottom], axis=0)
top = CB_G[0,] + np.array(pos_shift)
CB_G = np.append(CB_G, [top], axis=0)
#print(CB_G)

H2O = ss2.sdat
#delete particles from center
ss2.sdat.replicate_particles([3,3,4])
elem_num=len(ss2.sdat.particles['id'])
flag_2d = np.zeros((1,elem_num), dtype=int)
for c in CB_G:
    dist = (H2O.particles['pos'] - c)**2
    #同じ行のものを足す
    distarr=np.sum(dist, axis=1)
    flag = ((11000<distarr) & (distarr<12100))
    #bool -> int 
    flag = flag.astype(int)
    flag_insideCB = (distarr<10500)
    flag =  flag + (-1)*flag_insideCB.astype(int)  
    flag_2d = np.append(flag_2d, [flag], axis=0)
#print(flag_2d) 
flag = np.sum(flag_2d, axis=0) > 0
#print(flag)
ss2.sdat.trimming_particles(flag)

#trim over cell
flag = (0<H2O.particles['pos'][:,2]) & (H2O.particles['pos'][:,2]<cell_len[2])   
ss2.sdat.trimming_particles(flag)

#trim over cell
#Pt_flag = ss1.sdat.particles['type']==4
#Pt_pos = ss1.sdat.particles['pos'][Pt_flag]

Pt_list=[i for i in range(7,7+Pt_num)]

#this is for wrapped Pt
flag=((ss1.sdat.particles['mask']==10) & (ss1.sdat.particles['pos'][:,2]>cell_len[2]*2/3))
ss1.sdat.particles['pos'][flag]-=pos_shift
Pt_G = np.array([g_c_by_mask(ss1.sdat.particles, l) for l in Pt_list])

for pt_g in Pt_G:
    #trim H2O bulc close to Pt
    dist = (ss2.sdat.particles['pos'] - pt_g)**2
    distarr=np.sum(dist, axis=1)
    flag = distarr > 150
    ss2.sdat.trimming_particles(flag)
    #trimm terminal H or carbon
    dist = (ss1.sdat.particles['pos'] - pt_g)**2
    distarr=np.sum(dist, axis=1)
    dist_flag = distarr > 130
    type_flag = ss1.sdat.particles['type'] == 1
    flag = dist_flag & type_flag
    ss1.sdat.trimming_particles(flag)


#trim close Pt 
ss2.sdat.add_particles_property('mask', int, 1)
ss_out.sdat.concate_particles(ss2.sdat.particles)

print("number of water particles")
print(len(ss2.sdat.particles['id'])/3)
ss_out.output_file('test_inside_water.rd', 'input')
ss_out.output_file('show_dump', 'dumppos', ['id','type','pos','velo','mask'])

