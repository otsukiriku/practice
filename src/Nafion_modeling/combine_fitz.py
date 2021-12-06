import copy
import math
from MolCop import mmpystream as mmps
from MolCop.analysis import topology
import numpy as np
#construct object
CB_num=4
Pt_num=56

rmin = 110
rmax = 170
end_water = 124416

ss1 = mmps.Stream()
ss2 = mmps.Stream()

#import input name

ss1.import_file("trimmed_dump", 'dumppos')
#ss2.import_file("original_dump.pos.0", 'dumppos')
#ss2.import_file("original_dump.bond.0", 'dumpbond')
ss2.import_file("unwrapped_pos", 'dumppos')
#ss2.import_file("dump.bond.0", 'dumpbond')
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

def get_mol_gravity(atoms:mmps.Stream().sdat):
    G_list=[]
    m_list = list(set(atoms.particles['mol']))
    for ind, m in enumerate(m_list):
        flag = atoms.particles['mol'] == ind+1
        print(flag)
        gc = g_c(naf.particles, flag)
        #np.appendは遅いので，いったんpythonのリストにしてから再変換する
        gc = gc.tolist()
        G_list.append(gc[0])
    G_list = np.asarray(G_list)
    return G_list

def delete_mol_withG(atoms:mmps.Stream(), G_list=None):
#delete particles which Gravity is not in range rmin~rmax
    if G_list is None:
        print("you need to input gravity list (type:np.array)")
    else:
        #elem_num=len(atoms.particles['id'])
        m_list = list(set(atoms.particles['mol']))
        flag_1d = []
        flag_2d = []
        for c in CB_G:
            for g in G_list:
                dist = (g - c)**2
                #同じ行のものを足す
                distarr=np.sum(dist, axis=0)
                flag = (( (rmin**2) <distarr) & (distarr<(rmax**2)))
                #bool -> int 
                flag = flag.astype(int)
                flag_insideCB = (distarr<10500)
                flag_overcell = (g[2] > ss1.sdat.cell[2]) |(g[2] < ss1.sdat.dcell[2]) 
                flag =  flag + (-1)*flag_insideCB.astype(int) + (-1)*flag_overcell.astype(int)  
                flag = flag.tolist()
                flag_1d.append(flag)
            flag_2d.append(flag_1d)
        flag_2d = np.asarray(flag_2d)
        flag = np.sum(flag_2d, axis=0) > 0
        atoms.add_particles_property("flag", _dtype=bool, dim=1)
        #Nafionの重心だけを見たため，ナフィオンを構成する原子にflagをフィッティングする．
        for ind, m in enumerate(m_list):
            fl = atoms.particles['mol'] == ind+1
            atoms.particles['flag'][fl] = flag[ind]
        new_flag = atoms.particles['flag']
        atoms.trimming_particles(new_flag)


def delete_molecule(atoms:mmps.Stream(), d_mol=None):
    if d_mol is None:
        print("you need to input delete_molecule list")
    else:
        del_mol_flag=[]
        for i in del_mol:
            temp_flag = atoms.particles['mol'] == i
            temp_flag.tolist()
            del_mol_flag.append(temp_flag)
        del_mol_flag = np.asarray(del_mol_flag)
        del_mol_flag = np.sum(del_mol_flag.astype(int), axis = 0)
        flag = del_mol_flag == 0
        atoms.trimming_particles(flag)


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
"""
a=CB_num-1
bottom = CB_G[CB_num-1,] - np.array(pos_shift)
CB_G = np.append(CB_G, [bottom], axis=0)
top = CB_G[0,] + np.array(pos_shift)
CB_G = np.append(CB_G, [top], axis=0)
"""
#print(CB_G)
"""
flag = ss2.sdat.particles['id']<end_water+1
#print(flag)
H2O = copy.deepcopy(ss.sdat)
H2O.trimming_particles(flag)
H2O.replicate_particles([3,1,1])
"""
flag = ss2.sdat.particles['id']>end_water
naf = ss2.sdat
#naf.trimming_particles(flag)
#naf.replicate_particles([3,1,1])
#naf.create_connect_list(0.3)
#m_list = topology.create_molecule(naf)

#H2O_G=get_mol_gravity(H2O)
#delete_mol_withG(H2O, H2O_G)
naf_G=get_mol_gravity(naf)
delete_mol_withG(naf, naf_G)
print(naf_G)

#CBに近すぎるもの，入り込んでしまっているものを消す
flag_2d = []
for c in CB_G:
    dist = (naf.particles['pos'] - c)**2   
    dist = np.sum(dist, axis=1)
    flag_1d = dist < 10500
    flag_1d.tolist()
    flag_2d.append(flag_1d)
flag_2d = np.asarray(flag_2d)
flag = np.sum(flag_2d, axis=0) > 0
del_mol = naf.particles['mol'][flag]
del_mol = list(set(del_mol))
del_mol = np.asarray(del_mol)

delete_molecule(naf, d_mol=del_mol)


#Pt_flag = ss1.sdat.particles['type']==4
#Pt_pos = ss1.sdat.particles['pos'][Pt_flag]

Pt_list=[i for i in range(7,7+Pt_num)]

#this is for wrapped Pt
flag=((ss1.sdat.particles['mask']==10) & (ss1.sdat.particles['pos'][:,2]>cell_len[2]*2/3))
ss1.sdat.particles['pos'][flag]-=pos_shift

Pt_G = np.array([g_c_by_mask(ss1.sdat.particles, l) for l in Pt_list])

flag_2d=[]
np.zeros(0,)
naf.particles['flag']=np.zeros_like(naf.particles['id'], dtype=np.bool)
print(naf.particles['flag'])
for pt_g in Pt_G:
    #trim H2O bulc close to Pt
    naf_posx =  naf.particles['pos'][:,0]
    naf_posy =  naf.particles['pos'][:,1]
    naf_posz =  naf.particles['pos'][:,2]
    search_area_x = (pt_g[0]-50 < naf_posx) & (naf_posx < pt_g[0]+50)
    search_area_y = (pt_g[1]-50 < naf_posy) & (naf_posy < pt_g[1]+50)
    search_area_z = (pt_g[2]-50 < naf_posz) & (naf_posz < pt_g[2]+50)
    #print(len(naf_posx))
    search_area = search_area_x & search_area_y & search_area_z
    naf_pos = naf.particles['pos'][search_area]
    dist = (naf_pos - pt_g)**2
    distarr = np.sum(dist, axis=1)
    naf.particles['flag'][search_area] = distarr < 150
    flag_1d = naf.particles['flag']
    flag_1d.tolist()
    flag_2d.append(flag_1d)
#print(flag_2d)
flag_2d = np.asarray(flag_2d)
flag = np.sum(flag_2d, axis=0) > 0
del_mol = naf.particles['mol'][flag]
del_mol = list(set(del_mol))
del_mol = np.asarray(del_mol)

delete_molecule(naf, d_mol=del_mol)

#trim close Pt 
#ss2.sdat.add_particles_property('mask', int, 1)
ss_out.sdat.concate_particles(naf.particles)
ss_out.sdat.wrap_particles()

#print("number of water particles")
#print(len(ss2.sdat.particles['id'])/3)
#ss_out.output_file('test_inside_water.rd', 'input')
ss_out.output_file('show_dump', 'dumppos', ['id','type','pos','velo'])

