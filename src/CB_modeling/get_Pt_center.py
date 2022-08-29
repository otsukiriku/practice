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
from E_T.func import Pt_location,debug
from E_T.io import read_center as rc
import numpy as np
import sys
import copy
import math

#************************************
CB_num = 15
CB_radius = 75
Pt_num = 90
Pt_cluster_pnum = 309
Pt_mask_start = 16
#************************************
def g_c_by_mask(atoms : mmps.Stream().sdat.particles, mask):
    flag = atoms["mask"] == mask
    center = atoms["pos"][flag].mean(axis=0)
    return center

ss=mmps.Stream()
ss.import_file(sys.argv[1],'input')

#get each Ptcluster gravity center
Pt=ss.sdat
Pt_list=[i for i in range(Pt_mask_start,Pt_mask_start+Pt_num)]
Ptmask=Pt.particles['mask'][Pt.particles['type'] == 4]
Ptmask_2d=np.reshape(Ptmask, (Pt_num, Pt_cluster_pnum))

for idx, msk in enumerate(Ptmask_2d):
    #一列ずつPtmaskから切り取る．
    msk[:] = Pt_list[idx]
#print(Ptmask)

ss.sdat.particles['mask'][Pt.particles['type'] == 4] = Ptmask

Pt_list=np.array(Pt_list)
Pt_G = np.array([g_c_by_mask(Pt.particles, l) for l in Pt_list])


with open("Pt_center.txt","w"):
    print("Pt_center Ptnum:{}".format(len(Pt_G)))
    for pt_g,msk in zip(Pt_G,Pt_list):
        print("{} {} {} {}".format(msk,pt_g[0],pt_g[1],pt_g[2]))


