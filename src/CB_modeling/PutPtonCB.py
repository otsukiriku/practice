
from MolCop import mmpystream as mmps
from E_T.func import Pt_location
from E_T.io import read_center as rc
from E_T.io import import_rd as rd
from MolCop.analysis import neighbor
import random
import sys
import numpy as np
import math
sys.path.append("/nfshome12/rotsuki/molcop/src/")
import evaluate_structure.get_center as gc

#************************************
radii = 94
ptnum = 57
pt_clust_num = 309
cutoff = 60
Pt_mask_start = 11
#pr_tar = [10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,44,48,52]
pr_tar = [20,30,38,52]
#************************************

def construct_Pt(atoms:mmps.Stream().sdat,pt_pos):
    total_particle = pt_clust_num*ptnum

    if "id" not in atoms.particles:
        atoms.add_particles_property("id", _dtype=int)
    atoms.particles['id']=np.array([i+1 for i in range(total_particle)])

    if "pos" not in atoms.particles:
        atoms.add_particles_property("pos",dim=3)
    for _ in range(ptnum):
        atoms.particles['pos']=np.append(atoms.particles['pos'],pt_pos,axis=0)

    if "type" not in atoms.particles:
        atoms.add_particles_property("type", _dtype=int)
    atoms.particles['type']=np.array([4 for _ in range(total_particle)])

    if "mask" not in atoms.particles:
        atoms.add_particles_property("mask", _dtype=int)
    mask = []
    for i in range(Pt_mask_start, Pt_mask_start+ptnum):
        temp = [ i for _ in range(pt_clust_num)]
        mask.extend(temp)
    atoms.particles['mask']=np.array(mask)

def arrange_Pt_pos(atoms:mmps.Stream().sdat, Pt_mask, center, theta, phi):
    mask_flag = atoms.particles['mask'] == Pt_mask
    center_idx=random.randint(0,CB_num-1)
    theta[Pt_mask-Pt_mask_start]=random.uniform(0,360)
    phi[Pt_mask-Pt_mask_start]=random.uniform(0,360)
    ramda=random.uniform(0,360)
    theta = theta[Pt_mask-Pt_mask_start]
    phi = phi[Pt_mask-Pt_mask_start]

    x1 = atoms.particles['pos'][mask_flag][:,0]
    y1 = atoms.particles['pos'][mask_flag][:,1]*np.cos(phi)+atoms.particles['pos'][mask_flag][:,2]*np.sin(phi)
    z1 = -1*atoms.particles['pos'][mask_flag][:,1]*np.sin(phi)+atoms.particles['pos'][mask_flag][:,2]*np.cos(phi)
    x = center[center_idx][0] + x1*np.cos(theta) +  y1* np.sin(theta)
    y = center[center_idx][1] - x1*np.sin(theta) +  y1* np.cos(theta)
    z = center[center_idx][2] + z1
    shifted_pos = []
    [ shifted_pos.append([_x,_y,_z]) for _x, _y, _z in zip(x,y,z)]
    atoms.particles['pos'][mask_flag] = np.array(shifted_pos)

def check_dist_bet_Pt(Pt:mmps.Stream().sdat, Pt_mask):
    flag = True
    pt_flag = Pt.particles['mask']==Pt_mask
    center1 = gc.g_c(Pt.particles,pt_flag)
    #print(center1)
    for msk in range(Pt_mask-1,0,-1):
        fl = Pt.particles['mask']==Pt_mask
        center2=gc.g_c(Pt.particles,fl)
        #print(center2)
        dist = (center2 - center1)**2
        distarr = np.sqrt(np.sum(dist,axis=1))
        if distarr < cutoff :
            flag = False
            break
    return flag

def check_dist_bet_Pt2CB(Pt:mmps.Stream().sdat, Pt_mask, center):
    flag = True
    pt_flag = Pt.particles['mask']==Pt_mask
    center1 = gc.g_c(Pt.particles,pt_flag)
    #print(center1)
    dist = (center - center1)**2
    distarr = np.sqrt(np.sum(dist,axis=1))
    if (distarr < radii -2).any():
        flag = False
    return flag

if __name__ == "__main__":

    ss=mmps.Stream()
    Pt=mmps.Stream()
    ss.import_file(sys.argv[1],'dumppos')
    CB_num = len(set(ss.sdat.particles['mask']))
    CB_center = np.array(rc.r_c())
    Pt_pos,laxis=rd.import_ryudo(filename="Pt309.rd")

    Pt_pos = np.array(Pt_pos)
    laxis = np.array(laxis)

    Pt_pos = (Pt_pos - 0.5)*laxis*1e10
    Pt_pos[:,2]+=radii

    construct_Pt(Pt.sdat,Pt_pos)
    Pt.sdat.cell = ss.sdat.cell
    Pt.sdat.dcell = ss.sdat.dcell
    Pt.sdat.newcell = ss.sdat.newcell
    center_idx = np.array([i for i in range(CB_num)])
    theta = np.array([ 0 for i in range(ptnum)])
    phi = np.array([0 for i in range(ptnum)])

    end_flag = False
    #pr_flag = [False for _ in range(ptnum)]
    #for debug
    pr_flag = [True for _ in range(ptnum)]
    Pt_dist_flag = [True for _ in range(ptnum)]
    Pt2CB_dist_flag = [False for _ in range(ptnum)]
    while not end_flag:
        for idx,msk in enumerate(range(Pt_mask_start, Pt_mask_start+ptnum)):
            #print(msk)
            if pr_flag[idx] and Pt_dist_flag[idx] and Pt2CB_dist_flag[idx]:
                continue
            else:
                arrange_Pt_pos(Pt.sdat,msk,CB_center,theta,phi)
                """
                pr = Pt_location.Pt_loc_by_pos(Pt.sdat, msk, ptnum,Pt_mask_start=10)

                if idx < len(pr_tar) :
                    #print(pr_tar[idx])
                    if ( pr_tar[idx] -1 < pr) & ( pr < pr_tar[idx]+1):
                        pr_flag[idx] = True
                """
                #Pt_dist_flag[idx] = check_dist_bet_Pt(Pt.sdat,msk)
                Pt2CB_dist_flag[idx] = check_dist_bet_Pt2CB(Pt.sdat,msk,CB_center)
                if not ( pr_flag[idx] or Pt_dist_flag[idx] or Pt_dist_flag[idx]):
                    break
            #print("Pt_mask:{} probe radius:{}".format(msk,pr))

        #print("Pt_dist_flag")
        #print(Pt_dist_flag)
        #print("Pt2CB_dist_flag")
        #print(Pt2CB_dist_flag)
        if all(pr_flag) and all(Pt_dist_flag) and all(Pt2CB_dist_flag):
            end_flag = True

    print(Pt.sdat.particles['pos'])
    ss.sdat.concate_particles(Pt.sdat.particles)
    ss.sdat.total_particle += Pt.sdat.total_particle
    ss.output_file("show_dump","dumppos",["id","type","pos","mask"])

