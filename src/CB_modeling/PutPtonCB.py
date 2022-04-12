
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
CB_radii = 75
radii = CB_radii + 14
ptnum = 110
pt_clust_num = 309
cutoff = 60
Pt_mask_start = 11
#probe tar > 15
pr_tar = [15,15.5,16,16.5,17,17.5,18,18.5,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,56,58,60,62]
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

def check_dist_bet_Pt(atoms:mmps.Stream().sdat, Pt_mask):
    flag = True
    pt_flag = atoms.particles['mask']==Pt_mask
    center1 = gc.g_c(atoms.particles,pt_flag)
    #print(center1)
    Pt_location.correct_center_bymask(atoms,center1,Pt_mask)
    center1 = gc.g_c(atoms.particles,pt_flag)

    for msk in range(Pt_mask-1,Pt_mask_start-1,-1):
        #print(msk)
        fl = atoms.particles['mask']==msk
        center2=gc.g_c(atoms.particles,fl)
        dist = (center2 - center1)**2
        distarr = np.sum(dist,axis=1)
        #print(distarr)
        if distarr < cutoff**2 :
            flag = False
            break
    return flag

def check_dist_bet_Pt2CB(atoms:mmps.Stream().sdat, Pt_mask, center):
    flag = True
    pt_flag = atoms.particles['mask']==Pt_mask
    center1 = gc.g_c(atoms.particles,pt_flag)
    Pt_location.correct_center_bymask(atoms,center1,Pt_mask)
    center1 = gc.g_c(atoms.particles,pt_flag)
    #print(center1)
    dist = (center - center1)**2
    distarr = np.sum(dist,axis=1)
    closest_CB_idx = np.argmin(distarr)
    if (distarr < radii**2 - 4).any():
        flag = False
    Pt2CB_vector = center[closest_CB_idx] - center1
    atoms.particles['velo'][pt_flag] = Pt2CB_vector / 50000
    return flag

def check_Pt_G_posit_in_cell(atoms:mmps.Stream().sdat,Pt_mask):
    flag = True
    pt_flag = atoms.particles['mask']==Pt_mask
    center1 = gc.g_c(atoms.particles,pt_flag)
    Pt_location.correct_center_bymask(atoms,center1,Pt_mask)
    center1 = gc.g_c(atoms.particles,pt_flag)
    if (center1 > np.array(atoms.cell) - 10).any() or (center1 < np.array(atoms.dcell) + 10 ).any():
        flag = False
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
    Pt_G_flag = [False for _ in range(ptnum)]

    while not end_flag:
        for idx,msk in enumerate(range(Pt_mask_start, Pt_mask_start+ptnum)):
            #print(msk)
            if pr_flag[idx] and Pt_dist_flag[idx] and Pt2CB_dist_flag[idx] and Pt_G_flag[idx]:
                continue
            else:
                arrange_Pt_pos(Pt.sdat,msk,CB_center,theta,phi,Pt_pos)
                fl = Pt.sdat.particles['mask']==msk
                pr = Pt_location.Pt_loc_by_pos(Pt.sdat, msk, ptnum,Pt_mask_start=Pt_mask_start)

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
                Pt_G_flag[idx] = check_Pt_G_posit_in_cell(Pt.sdat,msk)
                if not Pt_G_flag[idx]:
                    break


        if all(pr_flag) and all(Pt_dist_flag) and all(Pt2CB_dist_flag) and all(Pt_G_flag):
            end_flag = True

    ss.sdat.concate_particles(Pt.sdat.particles)
    #ss.sdat.wrap_particles()
    ss.sdat.total_particle += Pt.sdat.total_particle
    Pt_list=[i for i in range(Pt_mask_start,Pt_mask_start+ptnum)]

    #Debug to show Probe
    pr_mask,pr_list,pr_pos = Pt_location.Pt_location(ss.sdat,CB_num=10,Pt_num=ptnum,Pt_mask_start=Pt_mask_start,prove_range=100)
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
    ss.sdat.wrap_particles()
    ss.sdat.set_elem_to_type(mmps.get_elem_to_type("para.rd"))
    ss.output_file("show_dump","dumppos",["id","type","pos",'velo',"mask"])

    Pt_flag=ss.sdat.particles['mask'] >= Pt_mask_start
    ss.sdat.particles['mask'][Pt_flag] = 11
    ss.sdat.thermofree_info = ['#thermofree 11']
    ss.output_file("newinput","input")

