from MolCop import mmpystream as mmps
import numpy as np
import sys
import math
import copy
from E_T.io import read_center as rc
ss = mmps.Stream()
flag_decomp = False
grid = np.zeros(3,dtype=int) #何回も利用する
totalgrid = 1
#*********************************************
#gridの幅
global grid_len
grid_len = 3
global CB_radii
CB_radii = 100
#*********************************************

def get_local_atoms(data : mmps.Stream().sdat):
    #gridの位置を設定
    dcmpx = list(range(0,grid[0], grid_len))
    dcmpy = list(range(0,grid[1], grid_len))
    dcmpz = list(range(0,grid[2], grid_len))
    #初期値0,x*y*zの3次元配列
    localatoms = np.zeros((grid[0],grid[1],grid[2]),dtype=int)
    atoms = data.particles
    for qx, qy, qz in atoms['pos']:
        i = binary_search(qx,dcmpx)*grid_len
        j = binary_search(qy,dcmpy)*grid_len
        k = binary_search(qz,dcmpz)*grid_len
        localatoms[i][j][k] += 1
    return localatoms

def inside_center(data:mmps.Stream().sdat,localatoms,CB_radius):
    CB_center = rc.r_c()
    cell = np.array(data.cell)
    dcell =  np.array(data.dcell)
    for center in CB_center:
        range_max = np.ceil((center + CB_radius)/grid_len)*grid_len
        range_min = np.floor((center - CB_radius)/grid_len)*grid_len
        over_flag = range_max >= cell
        under_flag = range_min < dcell
        if over_flag.any():
            range_max[over_flag] = cell[over_flag]
        if under_flag.any():
            range_min[under_flag] = dcell[under_flag]
        xrange= list(range(int(range_min[0]),int(range_max[0]), grid_len))
        yrange= list(range(int(range_min[1]),int(range_max[1]), grid_len))
        zrange= list(range(int(range_min[2]),int(range_max[2]), grid_len))
        for qz in zrange:
            for qy in yrange:
                for qx in xrange:
                    if ((qx-center[0])**2 + (qy-center[1])**2 + (qz-center[2])**2  < CB_radius**2 ):
                        localatoms[qx][qy][qz] += 1
    return localatoms

#２分探索(https://qiita.com/drken/items/97e37dd6143e33a64c8c)
def binary_search(key, arr):
    left = -1;
    right = len(arr)
    while(right - left > 1):
        mid = left + int((right - left) / 2)
        if(arr[mid] >= key): right = mid
        else: left = mid
    return right - 1

def set_grid(data,grid_len):
    global flag_decomp
    flag_decomp = data.flagdecomp
    global grid
    global totalgrid
    for i in range(3):
        grid[i] = math.ceil(data.cell[i])
        totalgrid *= math.ceil(grid[i]/grid_len)

def output_porosity(localatoms):
    filledgrid = np.count_nonzero(localatoms)
    porosity = (totalgrid-filledgrid)/totalgrid
    print("Totalgrid:{} * {} * {}  = {}".format(grid[0],grid[1],grid[2],totalgrid))
    print("Filledgrid:{}".format(filledgrid))
    print("Porosity:{}".format(porosity))

def debug_show(data:mmps.Stream().sdat,localatoms):
    debug = mmps.Stream()
    debug.sdat.cell = ss.sdat.cell
    debug.sdat.dcell = ss.sdat.dcell
    debug.sdat.newcell = ss.sdat.newcell
    new_pos = []
    dcmpx = list(range(0,grid[0], grid_len))
    dcmpy = list(range(0,grid[1], grid_len))
    dcmpz = list(range(0,grid[2], grid_len))
    for qz in dcmpz :
        for qy in dcmpy :
            for qx in dcmpx :
                if localatoms[qx][qy][qz] > 0:
                    new_pos.append([qx,qy,qz])
    new_pos = np.array(new_pos)
    total_particle = len(new_pos)
    debug.sdat.total_particle=total_particle

    if "id" not in debug.sdat.particles:
        debug.sdat.add_particles_property("id", _dtype=int)
    debug.sdat.particles['id']=np.array([i+1 for i in range(total_particle)])

    if "pos" not in debug.sdat.particles:
        debug.sdat.add_particles_property("pos",dim=3)
    debug.sdat.particles['pos']=new_pos

    if "type" not in debug.sdat.particles:
        debug.sdat.add_particles_property("type", _dtype=int)
    debug.sdat.particles['type']=np.array([10 for _ in range(total_particle)])
    
    debug.output_file("debug_show",'dumppos',['id','type','pos'])

if __name__ == '__main__':
    args = sys.argv
    inputfile = args[1]
    ss.import_file(inputfile, 'input')
    print("# Success : file import has completed !") #ファイル読み込み完了
    print("calculating now ...")
    set_grid(ss.sdat,grid_len)
    latoms = get_local_atoms(ss.sdat)
    #output_porosity(latoms)
    latoms = inside_center(ss.sdat,latoms, CB_radii)
    output_porosity(latoms)

    #debug
    debug_show(ss.sdat,latoms)

#############################################################
#遅い方法
    """
    for qx, qy, qz in atoms['pos']:
        flag = False
        for i in range(grid[0]):
            for j in range(grid[1]):
                for k in range(grid[2]):
                    minx = dcmpx[i]
                    maxx = dcmpx[i+1]
                    miny = dcmpy[j]
                    maxy = dcmpy[j+1]
                    minz = dcmpz[k]
                    maxz = dcmpz[k+1]
                    if minx<qx and qx<maxx and miny<qy and qy<maxy and minz<qz and qz<maxz:
                        localatoms[i][j][k] += 1;
                        flag = True
                        break
                if flag:
                    break
            if flag:
                break
    """
