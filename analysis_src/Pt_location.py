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
CB_num = 2
Pt_num = 20
Pt_cluster_pnum = 309

ss=mmps.Stream()
ss1=mmps.Stream()
#dim = 3 はxyzの要素数が3つあるので
ss1.sdat.add_particles_property('id',int,dim=1)
ss1.sdat.add_particles_property('pos',dim=3)
ss1.sdat.add_particles_property('type',int,dim=1)
ss.import_file('ptoncb.rd','input')
#ss.import_file('dump.bond.0','dumpbond')

def g_c_by_mask(atoms : mmps.Stream().sdat.particles, mask):
    flag = atoms["mask"] == mask
    center = atoms["pos"][flag].mean(axis=0)
    return center

#atom1 should be larger one
def get_between_dis(atom1, atom2):
#    for i in atom2
    dis = (atom1 - atom2)**2
    disarr = np.sum(dis,axis = 1)
    return disarr

# get center
dcell=ss.sdat.cell[2]
pos_shift = [0, 0, dcell]
CB_list=[i+1 for i in range(CB_num)]
CB_G = np.empty((0,3))
for mask in CB_list:
    #for wrap particles
    CB = ss.sdat.particles
    flag = ((CB['mask']==1) & (CB['pos'][:,2]>dcell*2/3))
    CB['pos'][flag]-=pos_shift
    flag = ((CB['mask']==CB_num) & (CB['pos'][:,2]<dcell*1/3))
    CB['pos'][flag]+=pos_shift

    flag = CB["mask"] == mask
    CB_g = [CB["pos"][flag].mean(axis=0)]
    CB_G = np.append(CB_G, CB_g, axis=0)
#print("CB_G")
#print(CB_G)

#get each Ptcluster gravity center
Pt=ss.sdat
Pt.create_connect_list(0.3)
#c_list = Pt.connect_list
#m_list = topology.create_molecule(Pt)
Pt_list=[i for i in range(7,7+Pt_num)]
#print(Pt_list)
Ptmask=Pt.particles['mask'][Pt.particles['type'] == 4]
Ptmask=np.reshape(Ptmask, (Pt_num, Pt_cluster_pnum))

for idx, msk in enumerate(Ptmask):
    #一列ずつPtmaskから切り取る．
    msk[:] = Pt_list[idx]
#print(Ptmask)

Ptmask = Ptmask.flatten()
ss.sdat.particles['mask'][Pt.particles['type'] == 4] = Ptmask

Pt_list=np.array(Pt_list)
print(Pt_list)
Pt_G = np.array([g_c_by_mask(Pt.particles, l) for l in Pt_list])
#print("Pt_G")
#print(Pt_G)

#get closest & 2nd closest dis from center
m_l1 = [i+1 for i in range(0,CB_num)]
mask_l = [m_l1 for j in range(Pt_num)]
Pt_CB_dis=np.array([get_between_dis(j, CB_G) for j in Pt_G])
#print("Pt_CB_dis")
#print(Pt_CB_dis)

#arguments of mask of closest CB from Pt
argsort = np.argsort(Pt_CB_dis, axis=1)

s_mask=np.empty((0,1), dtype=int)
s_dis=np.empty((0,1))
l_mask=np.empty((0,1),dtype=int)
l_dis=np.empty((0,1))
c_dis=np.empty((0,1))
for idx, arg in enumerate(argsort):
    #各Ptから最小のCBのmaskを取り出す．
    s_mask_idx = np.where(arg==0)[0][0]
    s_dis = np.append(s_dis,Pt_CB_dis[idx,s_mask_idx])
    s_mask = np.append(s_mask, s_mask_idx)
    #二番目に距離が短いmaskとかを取り出す．
    l_mask_idx = np.where(arg==1)[0][0]
    l_dis = np.append(l_dis,Pt_CB_dis[idx,l_mask_idx])
    l_mask = np.append(l_mask, l_mask_idx)
    #近い二つCBの中心間距離c_dis
    dis = ( CB_G[s_mask_idx] - CB_G[l_mask_idx])**2
    disarr = np.sum(dis,axis = 0)
    c_dis = np.append(c_dis, disarr)
    #Ptの重心から最近接，第二近接CBの中間の座標までの距離を出す
    CB_between_pos = ( CB_G[s_mask_idx] + CB_G[l_mask_idx] )/2
    Pt_to_CBbet = CB_between_pos - Pt_G

s_dis=np.sqrt(s_dis)
l_dis=np.sqrt(l_dis)
c_dis=np.sqrt(c_dis)
print("s_mask")
print(s_mask)
#print("s_dis")
#print(s_dis)
#print("l_dis")
#print(l_dis)
#print("c_dis")
#print(c_dis)

#prove radius
CB_to_Pt_vec=((Pt_G-CB_G[s_mask]))
pr = s_dis*c_dis**2/(-l_dis**2+s_dis**2+c_dis**2) - s_dis
correct = np.sum(s_dis)/Pt_num-100
pr = pr + correct
#print("prove radius")
#print(pr)

#get prove center(for test)
scalor=(s_dis+pr)/s_dis
scalor=np.reshape(scalor, (Pt_num,1))
pr_pos = CB_G[s_mask] + CB_to_Pt_vec*(scalor)
#print("pr_pos")
#print(pr_pos)
"""
print("cos theita")
print((-l_dis**2+s_dis**2+c_dis**2)/(2*s_dis*c_dis))
print("theita")
print(np.arccos((-l_dis**2+s_dis**2+c_dis**2)/(2*s_dis*c_dis))*180/3.14)
print("s_dis/(-l_dis**2+s_dis**2+c_dis**2)")
print(s_dis/(-l_dis**2+s_dis**2+c_dis**2))
"""
count = []
[count.append(i) for i in range(0,100,10)]
freq = []
mask = []
#Pt2CBbet = []
#for showing probe
pr_list = []
pr_pos_list = []
for co in count:
    pr_num=0
    for idx, r in enumerate(pr):
        if co <= r < co+10:
            pr_num+=1
            mask.append(Pt_list[idx])
            pr_list.append(pr[idx])
            pr_pos_list.append(pr_pos[idx])
            #Pt2CBbet.append(Pt_to_CBbet[idx])
            ss1.sdat.particles['id']=np.append(ss1.sdat.particles['id'],idx)
            ss1.sdat.particles['pos']=np.append(ss1.sdat.particles['pos'],[pr_pos[idx]],axis=0)
            ss1.sdat.particles['type']=np.append(ss1.sdat.particles['type'],Pt_list[idx])
    freq.append(pr_num)
#print(freq)
print("Ptmask_list")
print(mask)
print("pr_radius_list")
print(pr_list)
print("pr_pos_list")
[print(i) for i in pr_pos_list]
#print("Pt to CB between pos")
#[print(i) for i in Pt2CBbet]
print("prove_radius_range freq")
[print(f'{f}~{f+10} {c}', sep=' ', end="\n" ) for f, c in zip(count, freq)]

#print(ss.sdat.particles['pos'])


"""debug show probe radius
ss.sdat.concate_particles(ss1.sdat.particles)
print(ss.sdat.particles['id'])
#ss.output_file('test_input', 'input')
ss.output_file('probe_dump', 'dumppos',['id', 'type', 'pos'])
"""
