"""

ss.sdat <- input.rdを取得
center <- centerを取得
CB_G : 炭素担体の重心リスト
rc : 隣り合うマスク同士の中心間距離

Ptのmaskを一番近いセンターのmaskに合わせる
⇒mask は一番近い炭素担体のmaskが10の位, 次に近い炭素担体のmaskを1の位にする．
Pt_G : それぞれのPtの重心を出す．

rs : shorter distance from center
rs = sqrt(PtG^2 - center[mask])
rl : longer distance from center

mask 2~5
rl1 =sqrt(PtG^2 - center[mask+1])
rl2 =sqrt(PtG^2 - center[mask-1])
rl = bigger(rl1 or rl2)

rp : prove radius
rp = (rs * rc^2)/(-rl^2+rs^2+rc^2) -rs

"""


from MolCop import mmpystream as mmps
from MolCop.analysis import topology
from evaluate_structure import get_center as gc
import numpy as np
import sys
import copy
import math

ss=mmps.Stream()

ss.sdat.import_file('input.rd','input')

# get center
#dcell=ss.sdat.cell[2]
#pos_shift = [0, 0, -dcell]
#print(pos_shift)
for mask in range(1,7):
     = copy.deepcopy(ss.sdat.particles['mask']==mask)
    """
    # for wrapping particles
    flag = ((sdat_mask.particles['mask']==1) & (sdat_mask.particles['pos'][:,2]>dcell-30) & )
    sdat_mask.particles['pos'][flag]+=pos_shift
    """
    #get each CB gravity center
    CB.trimming_particles(sdat_mask.particles['type']==1)
    CB_pos = CB.particles['pos']
    CB_G = gc.g_c(CB_pos)
    #calc each CB distance in 2 dim list

CB_dis

#get each Ptcluster gravity center
Pt=copy.deepcopy(ss.sdat.particles['type']==4)
Pt.create_connect_list(0.3)
c_list = Pt.connect_list
m_list = topology.create_molecule(Pt)
for l in m_list:
    #get center
    flag = Pt.particles['mol'] == l
    Pt.trimming_particles[flag]
    Pt_G = gc.g_c(Pt.particles['pos'])
    #get closest & 2nd closest dis from center
    Pt.add_particles_property("close_dis")
    Pt.add_particles_property("close_mask")
    Pt.add_particles_property("2nd_close_dis")
    Pt.add_particles_property("2nd_close_mask")
    for i in CB_G:
        Pt_CB_pos = (Pt_G - i)**2
        Pt_CB_dis = np.sum(Pt_CB_pos, axis=1)
    # get closest CB's mask
    Pt.particles['close_mask'] = Pt_CB_dis.index(max(Pt_CB_dis))
    # pop closest Pt_CB_dis
    Pt.particles['close_dis'] = math.sqrt(Pt_CB_dis.pop(Pt.particles['close_mask']))
    Pt.particles['2nd_close_mask'] = Pt_CB_dis.index(max(Pt_CB_dis))
    Pt.particles['2nd_close_dis'] = math.sqrt(Pt_CB_dis.pop(Pt.particles['2nd_close_mask']))


    #print(mask - 1, end=' ')
    #[print(co, end=' ') for co in c_out]
    #print(end='\n')





