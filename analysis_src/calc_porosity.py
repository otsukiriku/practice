from MolCop import mmpystream as mmps
import numpy as np
import sys
import math
from E_T.io import read_center as rc
ss = mmps.Stream()
flag_decomp = False
grid = np.zeros(3,dtype=int) #何回も利用する
totalgrid = 1

def get_local_atoms(data : mmps.Stream().sdat):
    #gridの位置を設定
    dcmpx = [x for x in range(grid[0]+1)]
    dcmpy = [y for y in range(grid[1]+1)]
    dcmpz = [z for z in range(grid[2]+1)]

    #初期値0,x*y*zの3次元配列
    localatoms = np.zeros((grid[0],grid[1],grid[2]),dtype=int)
    atoms = data.particles
    #atomsの所属判定(二分探索)
    for qx, qy, qz in atoms['pos']:
        i = math.floor(qx)
        j = math.floor(qy)
        k = math.floor(qz)
        localatoms[i][j][k] += 1
    return localatoms

def inside_center(data:mmps.Stream().sdat,localatoms,CB_radius):
    CB_center = rc.r_c()
    cell = np.array(data.cell)
    dcell =  np.array(data.dcell)
    for center in CB_center:
        range_max = np.ceil(center + CB_radius)
        range_min = np.floor(center - CB_radius)
        over_flag = range_max >= cell
        under_flag = range_min < dcell
        if over_flag.any():
            range_max[over_flag] = cell[over_flag]
        if under_flag.any():
            range_min[under_flag] = dcell[under_flag]

        for qz in range(int(range_min[2]),int(range_max[2])):
            for qy in range(int(range_min[1]),int(range_max[1])):
                for qx in range(int(range_min[0]),int(range_max[0])):
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

def set_grid(data):
    global flag_decomp
    flag_decomp = data.flagdecomp
    global grid
    global totalgrid
    for i in range(3):
        grid[i] = math.ceil(data.cell[i])
        totalgrid *= math.ceil(grid[i])

def output_porosity(localatoms):
    filledgrid = np.count_nonzero(localatoms)
    porosity = (totalgrid-filledgrid)/totalgrid
    print("Totalgrid:{} * {} * {} = {}".format(grid[0],grid[1],grid[2],totalgrid))
    print("Filledgrid:{}".format(filledgrid))
    print("Porosity:{}".format(porosity))

#mainの時のみload_imbalanceを出力
if __name__ == '__main__':
    args = sys.argv
    inputfile = args[1]
    ss.import_file(inputfile, 'input')
    print("# Success : file import has completed !") #ファイル読み込み完了
    print("calculating now ...")
    set_grid(ss.sdat)
    latoms = get_local_atoms(ss.sdat)
    latoms = inside_center(ss.sdat,latoms, CB_radius=75)
    output_porosity(latoms)

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
