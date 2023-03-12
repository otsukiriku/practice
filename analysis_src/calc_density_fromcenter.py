from MolCop import mmpystream as mmps
import matplotlib.pyplot as plt
import numpy as np
import sys
import math
import copy
from E_T.io import read_center as rc
from E_T.func import debug
#*********************************************
CB_radii = 75
CB_num= 15
line_color=["blue","red", "green", "orange","b"]
#*********************************************

def density_from_center(atoms:mmps.Stream().sdat, tar_type, center):
    count = np.zeros((1,1000),dtype=int)[0]
    target = atoms.particles["pos"][atoms.particles['type']==tar_type]
    print(len(target))
    for c in center:
        distarr = (target - c)**2
        distarr = np.round(np.sqrt(np.sum(distarr,axis=1)))
        distarr = distarr.astype(int)
    print(distarr.shape)
    for d in distarr:
            count[d]+=1
    return count

def create_graph(xarray,yarray,data_title,figtitle):
    """x,y should be 2 dim array"""
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)
    ax.set_ylim(0, 0.0005)
    ax.set_xlabel("r(Å)", size = 16,)
    ax.set_ylabel('Density', size = 16,)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    #plt.figure(figsize=(4,3))
    for idx,(x,y,d_t) in enumerate(zip(xarray,yarray,data_title)):
        ax.plot(x,y,color=line_color[idx])
        ax.text(0.95,0.75+0.075*idx, line_color[idx]+":"+d_t, transform=ax.transAxes, horizontalalignment="right",size=20)
        fig.savefig( figtitle + ".png",dpi=250)

def g_c(atoms : mmps.Stream().sdat, flag=None):
    if flag is None:
        center = np.array(atoms.particles['pos'].mean(axis=0))
    else:
        center = np.array(atoms.particles['pos'][flag].mean(axis=0))
    return center

def correct_center(atoms : mmps.Stream().sdat, center:np.array):
    half_cell=CB_radii +100
    for idx,cent in enumerate(center):
        flag=ss.sdat.particles['mask']==idx+1
        over_flag = (atoms.particles['pos'][flag]-cent) > half_cell
        under_flag = (atoms.particles['pos'][flag]-cent) < -half_cell
        atoms.particles['pos'][flag] -= over_flag * atoms.newcell
        atoms.particles['pos'][flag] += under_flag * atoms.newcell


if __name__ == '__main__':
    #data_num = len(sys.argv) - 1
    sys.argv.pop(0)
    data_title = sys.argv
    print(data_title)
    dist_array = []
    density_array = []
    corrected_center = []
    for d_t in data_title:
        ss = mmps.Stream()
        ss.import_file(d_t, 'dumppos')
        CB_center=[]
        for mask in range(1,CB_num+1):
            fl=ss.sdat.particles['mask']==mask
            c_out = g_c(ss.sdat,fl)
            CB_center.append(c_out)
        CB_center=np.array(CB_center)
        correct_center(ss.sdat,CB_center)
        CB_center=[]
        for mask in range(1,CB_num+1):
            fl=ss.sdat.particles['mask']==mask
            c_out = g_c(ss.sdat,fl)
            CB_center.append(c_out)
        CB_center=np.array(CB_center)
        debug.visualize_by_add_atom(ss.sdat,CB_center)
        ss.output_file("gravity"+d_t, "dumppos",["id","type","pos"])

        count = density_from_center(ss.sdat, 1, CB_center)
        particle_num=len(CB_center)
        temp_dist = []
        temp_density = []
        for i,c in enumerate(count):
            if not c == 0:
                density = c/(4.0/3.0*math.pi*((i+1)**3-i**3))/particle_num
                temp_dist.append(i)
                temp_density.append(density)
                #print("{} {}".format(i,density))
            if i == 200:
                break
        dist_array.append(temp_dist)
        density_array.append(temp_density)
    create_graph(dist_array,density_array,data_title=data_title,figtitle="STEP0vs500000vs1000000_densityfromcenter")
