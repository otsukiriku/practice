"""

ss.sdat <- input.rdを取得
center <- centerを取得
CB_G : 炭素担体の重心リスト(二次元)[[x,y,z], ....]

Ptのmaskを一番近いセンターのmaskに合わせる
⇒mask は一番近い炭素担体のmaskが10の位, 次に近い炭素担体のmaskを1の位にする．
Pt_G : それぞれのPtの重心を出す．

s_dis(rs) : shorter distance from center list
l_dis(rl) : longer distance from center list
c_dis(rc) : 隣り合うマスク同士の中心間距離(1次元)

mask 2~5
rl1 =sqrt(PtG^2 - center[mask+1])
rl2 =sqrt(PtG^2 - center[mask-1])
rl = bigger(rl1 or rl2)

rp : prove radius
rp = (rs * rc^2)/(-rl^2+rs^2+rc^2) -rs

"""


from MolCop import mmpystream as mmps
from MolCop.analysis import topology
#from evalueate_structure import get_center as gc
import numpy as np
import sys
import copy
import math

ss=mmps.Stream()

ss.import_file('ptoncb.rd','input')
ss.import_file('dump.bond.20000','dumpbond')

def g_c(atoms : mmps.Stream().sdat.particles):
    x, y, z = np.hsplit(atoms['pos'],3)
    center = np.array([x.mean(), y.mean(), z.mean()])
    return center

def get_between_dis(atom1, atom2):
    dis = (atom1 - atom2)**2
    disarr = np.sum(dis)
    return disarr
# get center
dcell=ss.sdat.cell[2]
pos_shift = [0, 0, -dcell]
#print(pos_shift)
CB_G = np.empty((0,3))
for mask in range(1,7):
    CB = copy.deepcopy(ss.sdat)
    CB.trimming_particles(CB.particles['mask']==mask)
    #print(f'CB position {CB.particles["pos"]}')
    """
    # for wrapping particles
    flag = ((sdat_mask.particles['mask']==1) & (sdat_mask.particles['pos'][:,2]>dcell-30) & )
    sdat_mask.particles['pos'][flag]+=pos_shift
    """
    #get each CB gravity center
    flag = (CB.particles['type']==1)
    CB.trimming_particles(flag)
    flag = ((CB.particles['mask']==1) & (CB.particles['pos'][:,2]>750))
    CB.particles['pos'][flag]+=pos_shift

    CB_pos = CB.particles
    #print(f'center{CB_pos}')
    CB_g = g_c(CB_pos)
    #print(CB_g)
    CB_G = np.append(CB_G, [CB_g], axis=0)

    #calc each CB distance in 2 dim list


for x in CB_G:
    for y in CB_G:
        CB_Gpos = get_between_dis(x, y)
        #print(f'CB gravity center position{CB_Gpos}')


#get each Ptcluster gravity center
Pt=copy.deepcopy(ss.sdat)
print(Pt.bondorder_list)
flag = Pt.particles['type']==4
print(flag)
Pt.trimming_particles(flag)
Pt.create_connect_list(0.3)
print(Pt.connect_list)

c_list = Pt.connect_list
m_list = topology.create_molecule(Pt)



'''
ss.sdat.add_particles_property("molnum", _dtype = int, dim =1)
for l in m_list:
    #get center
    num_mol = len(l)
    for ind in l:
        ss.sdat.particles["molnum"][ind]=num_mol
    print(ss.sdat.particles["molnum"])
    flag = Pt.particles['mol'] == l
    Pt.trimming_particles[flag]
    Pt_G = g_c(Pt.particles)
    #get closest & 2nd closest dis from center
    """
    Pt.add_particles_property("close_dis")
    Pt.add_particles_property("close_mask")
    Pt.add_particles_property("2nd_close_dis")
    Pt.add_particles_property("2nd_close_mask")
    """
    for i in CB_G:
        Pt_CB_Gpos = (Pt_G - i)**2
        Pt_CB_dis = np.sum(Pt_CB_Gpos, axis=1)
    # get closest CB's mask
    s_mask = Pt_CB_dis.index(max(Pt_CB_dis))
    # pop closest Pt_CB_dis
    s_dis = math.sqrt(Pt_CB_dis.pop(s_mask))
    l_mask = Pt_CB_dis.index(max(Pt_CB_dis))
    l_dis = math.sqrt(Pt_CB_dis.pop(l_mask))
    c_dis = CB_dis[s_mask][l_mask]


#prove radius
pr = (s_dis * c_dis**2)/(-l_dis**2+s_dis**2+c_dis**2) - s_dis

count = 0
for count in range(0,50,2):
    if (count < pr)&(pr < count+2):
        freq[count]+=1


[print(count, freq[count], sep=' ', end="\n" ) for count in range(0,50,2)]
file = open("./Pt_location.txt", encofing = 'UTF-8')
print(count, freq[count], sep=' ', end="\n" )
file.close

#print(mask - 1, end=' ')
    #[print(co, end=' ') for co in c_out]
    #print(end='\n')




'''
