
from MolCop import mmpystream as mmps
from E_T.func import Pt_location,get_center,debug
from E_T.io import read_center as rc
from E_T.io import import_rd as rd
from MolCop.analysis import neighbor
import random
import sys
import numpy as np
import math
import copy
sys.path.append("/nfshome12/rotsuki/molcop/src/")
import evaluate_structure.get_center as gc

#************************************
radii = 94
ptnum = 57
pt_clust_num = 309
cutoff = 60
Pt_mask_start = 11
pr_tar = [10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,44,48,52]
#pr_tar = [20,30,38,52]
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

    if "velo" not in atoms.particles:
        atoms.add_particles_property("velo",dim=3)
    for _ in range(ptnum):
        atoms.particles['velo']=np.append(atoms.particles['velo'],pt_pos,axis=0)

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

def arrange_Pt_pos(atoms:mmps.Stream().sdat, Pt_mask, center, theta, phi, Init_Ptpos):
    mask_flag = atoms.particles['mask'] == Pt_mask
    center_idx=random.randint(0,CB_num-1)
    theta[Pt_mask-Pt_mask_start]=random.uniform(0,360)
    phi[Pt_mask-Pt_mask_start]=random.uniform(0,360)
    ramda=random.uniform(0,360)
    theta = theta[Pt_mask-Pt_mask_start]
    phi = phi[Pt_mask-Pt_mask_start]

    x1 = Init_Ptpos[:,0]
    y1 = Init_Ptpos[:,1]*np.cos(phi)+Init_Ptpos[:,2]*np.sin(phi)
    z1 = -1*Init_Ptpos[:,1]*np.sin(phi)+Init_Ptpos[:,2]*np.cos(phi)
    x = center[center_idx][0] + x1*np.cos(theta) +  y1* np.sin(theta)
    y = center[center_idx][1] - x1*np.sin(theta) +  y1* np.cos(theta)
    z = center[center_idx][2] + z1
    shifted_pos = []
    [ shifted_pos.append([_x,_y,_z]) for _x, _y, _z in zip(x,y,z)]
    atoms.particles['pos'][mask_flag] = np.array(shifted_pos)
    #wrap at boundary
    over_cell = atoms.particles['pos'][mask_flag] >= np.array(atoms.cell)
    under_cell = atoms.particles['pos'][mask_flag] < np.array(atoms.dcell)
    atoms.particles['pos'][mask_flag] -= over_cell * np.array(atoms.newcell)
    atoms.particles['pos'][mask_flag] += under_cell * np.array(atoms.newcell)

def replicate_center_xyz(center:np.array,cell):
    center_org = copy.deepcopy(center)
    cell= np.array(cell)
    wrap_list=np.array([[1,0,0],[0,1,0],[0,0,1]])
    for lst in wrap_list:
        p_shift=center_org+cell*lst
        m_shift=center_org-cell*lst
        center=np.append(center,p_shift,axis=0)
        center=np.append(center,m_shift,axis=0)
    return center

""" mistake
    cell = np.array(cell)
    org_center = copy.deepcopy(center)
    rep_center = copy.deepcopy(center)
    for idx,c in enumerate(cell):
        plus_center=org_center
        minus_center=org_center
        shift = np.array([0,0,0])
        shift[idx] = c
        np.add(plus_center,shift)
        np.add(minus_center,shift)
        print(center)
        #minus_center -= shift
        rep_center = np.append(rep_center,plus_center,axis=0)
        rep_center = np.append(rep_center,minus_center,axis=0)
    return rep_center
"""

def check_dist_bet_Pt(Pt:mmps.Stream().sdat, Pt_mask):
    flag = True
    pt_flag = Pt.particles['mask']==Pt_mask
    center1 = gc.g_c(Pt.particles,pt_flag)
    #print(center1)
    Pt_location.correct_center_bymask(Pt,center1,Pt_mask)
    center1 = gc.g_c(Pt.particles,pt_flag)

    for msk in range(Pt_mask-1,Pt_mask_start,-1):
        #print(msk)
        fl = Pt.particles['mask']==msk
        center2=gc.g_c(Pt.particles,fl)
        dist = (center2 - center1)**2
        distarr = np.sqrt(np.sum(dist,axis=1))
        #print(distarr)
        if distarr < cutoff :
            flag = False
            break
    return flag

def check_dist_bet_Pt2CB(Pt:mmps.Stream().sdat, Pt_mask, center):
    flag = True
    pt_flag = Pt.particles['mask']==Pt_mask
    center1 = gc.g_c(Pt.particles,pt_flag)
    Pt_location.correct_center_bymask(Pt,center1,Pt_mask)
    center1 = gc.g_c(Pt.particles,pt_flag)
    #print(center1)
    dist = (center - center1)**2
    distarr = np.sqrt(np.sum(dist,axis=1))
    closest_CB_idx = np.argmin(distarr)
    if (distarr < radii -2).any():
        flag = False
    Pt2CB_vector = center[closest_CB_idx] - center1
    Pt.particles['velo'][pt_flag] = Pt2CB_vector / 5000
    return flag

if __name__ == "__main__":

    ss=mmps.Stream()
    Pt=mmps.Stream()
    ss.import_file(sys.argv[1],'dumppos')
    CB_num = len(set(ss.sdat.particles['mask']))
    CB_center = np.array(rc.r_c())
    ss_cell = ss.sdat.cell
    rep_CB_center = replicate_center_xyz(CB_center,ss_cell)
    Pt_pos,laxis=rd.import_ryudo(filename="Pt309.rd")

    #Read Pt position
    Pt_pos = np.array(Pt_pos)
    laxis = np.array(laxis)
    Pt_pos = (Pt_pos - 0.5)*laxis*1e10
    Pt_pos[:,2]+=radii

    #Preparation
    construct_Pt(Pt.sdat,Pt_pos)
    Pt.sdat.cell = ss.sdat.cell
    Pt.sdat.dcell = ss.sdat.dcell
    Pt.sdat.newcell = ss.sdat.newcell
    center_idx = np.array([i for i in range(CB_num)])
    theta = np.array([ 0 for _ in range(ptnum)])
    phi = np.array([ 0 for _ in range(ptnum)])
    end_flag = False
    pr_flag = [False for _ in range(ptnum)]
    Pt_dist_flag = [False for _ in range(ptnum)]
    Pt2CB_dist_flag = [False for _ in range(ptnum)]

    while not end_flag:
        for idx,msk in enumerate(range(Pt_mask_start, Pt_mask_start+ptnum)):
            #print(msk)
            if pr_flag[idx] and Pt_dist_flag[idx] and Pt2CB_dist_flag[idx]:
                continue
            else:
                arrange_Pt_pos(Pt.sdat,msk,CB_center,theta,phi,Pt_pos)
                fl = Pt.sdat.particles['mask']==msk
                pr = Pt_location.Pt_loc_by_pos(Pt.sdat, msk, ptnum,Pt_mask_start=Pt_mask_start)
                #pr = Pt_location.Pt_location(Pt.sdat,nonzero_option=False, CB_num=10, Pt_num=1,Pt_mask_start=msk)

                if idx < len(pr_tar) :
                    #print(pr_tar[idx])
                    if ( pr_tar[idx] -1 < pr) & ( pr < pr_tar[idx]+1):
                        pr_flag[idx] = True
                        print("mask:{}, probe:{}".format(msk,pr))
                    else:
                        break
                else:
                    if pr > max(pr_tar):
                        pr_flag[idx] = True
                    else:
                        break

                if idx == 0:
                    Pt_dist_flag[idx] = True
                else:
                    Pt_dist_flag[idx] = check_dist_bet_Pt(Pt.sdat,msk)
                    if not Pt_dist_flag[idx]:
                        break

                Pt2CB_dist_flag[idx] = check_dist_bet_Pt2CB(Pt.sdat,msk,rep_CB_center)
                if not Pt2CB_dist_flag[idx]:
                    break


        if all(pr_flag) and all(Pt_dist_flag) and all(Pt2CB_dist_flag):
            end_flag = True

    ss.sdat.concate_particles(Pt.sdat.particles)
    #ss.sdat.wrap_particles()
    ss.sdat.total_particle += Pt.sdat.total_particle
    Pt_list=[i for i in range(Pt_mask_start,Pt_mask_start+ptnum)]

    #Debug to show Probe
    pr_mask,pr_list,pr_pos = Pt_location.Pt_location(ss.sdat,CB_num=10,Pt_num=ptnum,Pt_mask_start=Pt_mask_start)
    [print("Ptmask:{}  Probe:{}".format(msk,pr)) for msk,pr in zip(pr_mask,pr_list)]
    pr_type = []
    #for msk in range(Pt_mask_start, Pt_mask_start+ptnum):
    [ pr_type.append(msk) for msk in range(Pt_mask_start,Pt_mask_start+ptnum) ]
    #debug.visualize_by_add_atom(ss.sdat, tar_pos=pr_pos, tar_type=np.array(pr_type))
    #print(ss.sdat.particles['id'].shape)
    #print(ss.sdat.particles['type'].shape)
    #print(ss.sdat.particles['pos'].shape)
    #print(ss.sdat.particles['velo'].shape)
    #print(ss.sdat.particles['mask'].shape)
    """
    #Debug to show Pt_G
    Pt_G = []
    for l in Pt_list:
        fl = Pt.sdat.particles['mask'] == l
        Pt_G.append(gc.g_c(Pt.sdat.particles,fl)[0])
    Pt_G = np.array(Pt_G)
    debug.visualize_by_add_atom(ss.sdat,Pt_G)
    """

    #debug All Pt_mask will be
    Pt_flag=ss.sdat.particles['mask'] >= Pt_mask_start
    ss.sdat.particles['mask'][Pt_flag] = 4

    ss.output_file("newinput","input")
    ss.output_file("show_dump","dumppos",["id","type","pos",'velo',"mask"])

