
from sympy import per
from E_T.io import read_center as rc
from glob import glob
from random import random, randrange
from MolCop import mmpystream as mmps
from MolCop.analysis import topology
import numpy as np
import sys
import math
import copy
sys.path.append("/nfshome12/rotsuki/molcop/src/")
ss = mmps.Stream()
flag_decomp = False
grid = np.zeros(3, dtype=int)  # 何回も利用する
totalgrid = 1
# *********************************************
# gridの幅
global grid_len
grid_len = 2
global CB_radii
CB_radii = 78
global O2_molecule_num
O2_molecule_num = 5900
global O2_outside_limit_from_center
O2_outside_limit_from_center = 130
# *********************************************


def get_local_atoms(data: mmps.Stream().sdat):
    # gridの位置を設定
    dcmpx = list(range(0, grid[0], grid_len))
    dcmpy = list(range(0, grid[1], grid_len))
    dcmpz = list(range(0, grid[2], grid_len))
    # 初期値0,x*y*zの3次元配列
    localatoms = np.ones((grid[0]+1, grid[1]+1, grid[2]+1), dtype=int)
    atoms = data.particles
    for qx, qy, qz in atoms['pos']:
        i = binary_search(qx, dcmpx)*grid_len
        j = binary_search(qy, dcmpy)*grid_len
        k = binary_search(qz, dcmpz)*grid_len
        localatoms[i][j][k] += 1
    return localatoms


def fill_o2_outside_limit_from_center(data: mmps.Stream().sdat, localatoms):
    CB_center = rc.r_c()
    cell = np.array(data.cell)
    dcell = np.array(data.dcell)
    for center in CB_center:
        range_max = np.ceil(
            (center + O2_outside_limit_from_center)/grid_len)*grid_len
        range_min = np.floor(
            (center - O2_outside_limit_from_center)/grid_len)*grid_len
        over_flag = range_max >= cell
        under_flag = range_min < dcell
        if over_flag.any():
            range_max[over_flag] = cell[over_flag]
        if under_flag.any():
            range_min[under_flag] = dcell[under_flag]
        xrange = list(range(int(range_min[0]), int(range_max[0]), grid_len))
        yrange = list(range(int(range_min[1]), int(range_max[1]), grid_len))
        zrange = list(range(int(range_min[2]), int(range_max[2]), grid_len))
        for qz in zrange:
            for qy in yrange:
                for qx in xrange:
                    if ((qx-center[0])**2 + (qy-center[1])**2 + (qz-center[2])**2 < O2_outside_limit_from_center**2):
                        localatoms[qx][qy][qz] -= 1
    return localatoms


def fill_inside_center(data: mmps.Stream().sdat, localatoms, CB_radius):
    CB_center = rc.r_c()
    cell = np.array(data.cell)
    dcell = np.array(data.dcell)
    for center in CB_center:
        range_max = np.ceil((center + CB_radius)/grid_len)*grid_len
        range_min = np.floor((center - CB_radius)/grid_len)*grid_len
        over_flag = range_max >= cell
        under_flag = range_min < dcell
        if over_flag.any():
            range_max[over_flag] = cell[over_flag]
        if under_flag.any():
            range_min[under_flag] = dcell[under_flag]
        xrange = list(range(int(range_min[0]), int(range_max[0]), grid_len))
        yrange = list(range(int(range_min[1]), int(range_max[1]), grid_len))
        zrange = list(range(int(range_min[2]), int(range_max[2]), grid_len))
        for qz in zrange:
            for qy in yrange:
                for qx in xrange:
                    if ((qx-center[0])**2 + (qy-center[1])**2 + (qz-center[2])**2 < CB_radius**2):
                        localatoms[qx][qy][qz] += 3
    return localatoms


# ２分探索(https://qiita.com/drken/items/97e37dd6143e33a64c8c)


def binary_search(key, arr):
    left = -1
    right = len(arr)
    while(right - left > 1):
        mid = left + int((right - left) / 2)
        if(arr[mid] >= key):
            right = mid
        else:
            left = mid
    return right - 1


def set_grid(data, grid_len):
    global flag_decomp
    flag_decomp = data.flagdecomp
    global grid
    global totalgrid
    for i in range(3):
        grid[i] = math.ceil(data.cell[i])
        totalgrid *= math.ceil(grid[i]/grid_len)


def put_O2(data: mmps.Stream().sdat, localatoms):
    ss1 = mmps.Stream()
    ss1.import_file("O2.dump", "dumppos")
    ss1.sdat.particles["pos"] -= ss1.sdat.particles["pos"][0]
    ss_tmp = copy.deepcopy(ss1)
    O2_data = mmps.Stream().sdat
    counter = O2_molecule_num
    while(counter > 0):
        rand_x = randrange(0, grid[0])
        rand_y = randrange(0, grid[1])
        rand_z = randrange(0, grid[2])
        if localatoms[rand_x][rand_y][rand_z] >= 1:
            continue
        else:
            localatoms[rand_x][rand_y][rand_z] += 1
            axis = np.array([randrange(0, 2), randrange(0, 2), 1], dtype=float)
            angle = float(randrange(0, 180))
            rotted_pos = Rotate_particles(ss1.sdat.particles, axis, angle)
            ss_tmp.sdat.particles["pos"] = rotted_pos + np.array(
                [rand_x+grid_len/2, rand_y + grid_len/2, rand_z+grid_len/2])
            O2_data.concate_particles(ss_tmp.sdat.particles)
            counter -= 1
            # if counter % 100 == 0:
            #     print("leftParticleNum:{:3d}".format(counter))
    return O2_data


def output_porosity(data: mmps.Stream().sdat, localatoms):
    filledgrid = np.count_nonzero(localatoms)
    porosity = (totalgrid-filledgrid)/totalgrid
    print(
        "Totalgrid:{} * {} * {}  = {}".format(grid[0], grid[1], grid[2], totalgrid))
    print("Filledgrid:{}".format(filledgrid))
    print("Porosity:{}".format(porosity))


def debug_show(data: mmps.Stream().sdat, localatoms):
    debug = mmps.Stream()
    debug.sdat.cell = ss.sdat.cell
    debug.sdat.dcell = ss.sdat.dcell
    debug.sdat.newcell = ss.sdat.newcell
    new_pos = []
    dcmpx = list(range(0, grid[0], grid_len))
    dcmpy = list(range(0, grid[1], grid_len))
    dcmpz = list(range(0, grid[2], grid_len))
    for qz in dcmpz:
        for qy in dcmpy:
            for qx in dcmpx:
                if localatoms[qx][qy][qz] > 0:
                    new_pos.append([qx, qy, qz])
    new_pos = np.array(new_pos)
    total_particle = len(new_pos)
    debug.sdat.total_particle = total_particle

    if "id" not in debug.sdat.particles:
        debug.sdat.add_particles_property("id", _dtype=int)
    debug.sdat.particles['id'] = np.array([i+1 for i in range(total_particle)])

    if "pos" not in debug.sdat.particles:
        debug.sdat.add_particles_property("pos", dim=3)
    debug.sdat.particles['pos'] = new_pos

    if "type" not in debug.sdat.particles:
        debug.sdat.add_particles_property("type", _dtype=int)
    debug.sdat.particles['type'] = np.array(
        [10 for _ in range(total_particle)])
    debug.output_file("debug_show", 'dumppos', ['id', 'type', 'pos'])


def Rodrigues(n, th):
    norm = math.sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2])
    if norm == 0.0:
        print("Error : invalid axis : ", n)
        sys.exit()
    n /= norm
    n_1d = n.reshape(1, -1)
    nT = n.reshape(-1, 1)
    R0 = math.cos(math.radians(th)) * np.eye(3)
    R1 = (1 - math.cos(math.radians(th))) * (nT @ n_1d)
    R2 = math.sin(math.radians(th)) * np.array([[0.0, -n[2], n[1]],
                                                [n[2], 0.0, -n[0]],
                                                [-n[1], n[0], 0.0]], dtype=float)
    R = R0 + R1 + R2
    return R


def Rotate_particles(atoms: mmps.Stream().sdat.particles, nvec, theta):
    R_n = Rodrigues(nvec, theta)
    center = atoms["pos"].mean(axis=0)
    atoms['pos'] -= center
    atoms['pos'] = atoms['pos'] @ R_n
    atoms['pos'] += center
    return atoms["pos"]


def separate_NafionData(data: mmps.Stream().sdat):
    Pt_idx = np.min(data.particles["id"][data.particles["type"] == 4])
    start_idx_ofNafion = np.min(
        data.particles["id"][data.particles["type"] == 1] > Pt_idx)
    flag = data.particles["id"] >= start_idx_ofNafion
    nafion_data = mmps.Stream().sdat
    nafion_data = copy.deepcopy(data)
    nafion_data.trimming_particles(flag)
    reversed_flag = np.logical_not(flag)
    data.trimming_particles(reversed_flag)
    return nafion_data


def create_molnum(data: mmps.Stream().sdat):
    data.create_connect_list(0.3)
    connect_list = data.connect_list
    m_list = topology.create_molecule(data)

    data.add_particles_property("molnum", _dtype=int, dim=1)
    for l in m_list:
        num_mol = len(l)
        for ind in l:
            data.particles["molnum"][ind] = num_mol


def delete_H_of_Nafion(data: mmps.Stream().sdat, percentage=80):
    start_idx_ofNafion = np.min(
        data.particles["id"][data.particles["type"] == 6])-28
    flag_Nafion = data.particles["id"] >= start_idx_ofNafion
    # print(start_idx_ofNafion)
    # print(ss.sdat.particles["id"][flag_Nafion])

    flag_hydrogen_from_Nafion = flag_Nafion & (data.particles["type"] == 2)
    random_boolean = np.array(
        [randrange(0, 100) <= percentage for i in range(ss.sdat.total_particle)])
    flag_delete = flag_hydrogen_from_Nafion & random_boolean

    num_of_deleted_hydrogen = np.sum(
        flag_hydrogen_from_Nafion) - np.sum(flag_delete)

    # print("全てのHの数")
    # print(np.sum(flag_hydrogen_from_Nafion))
    # print("置換するHの数")
    # print(np.sum(flag_delete))

    flag_delete = np.logical_not(flag_delete)
    data.trimming_particles(flag_delete)
    return num_of_deleted_hydrogen


def change_H2O_to_H3O(data: mmps.Stream().sdat, num_of_deleted_hydrogen):
    # type_lst = data.particles["type"][data.particles["molnum"] == 3]
    # type_list_with_molnum_is_3 = list(set(data.particles["mol"][data.particles["molnum"] == 3])
    molnum3_flag = data.particles["molnum"] == 3
    mol_lst = list(data.particles["mol"][molnum3_flag])
    mol_set = list(set(data.particles["mol"][molnum3_flag]))
    id_lst = list(data.particles["id"][molnum3_flag])
    type_lst = list(data.particles["type"][molnum3_flag])
    pos_lst = list(data.particles["pos"][molnum3_flag])

    dict_mol_id = dict()
    dict_mol_type = dict()
    dict_mol_pos = dict()
    for m, i, t, p in zip(mol_lst, id_lst, type_lst, pos_lst):
        if m not in dict_mol_type:
            dict_mol_id[m] = [i]
            dict_mol_type[m] = [t]
            dict_mol_pos[m] = [p]
        else:
            dict_mol_id[m].append(i)
            dict_mol_type[m].append(t)
            dict_mol_pos[m].append(p)

    # for k,v in dict_mol_type.items():
    #     print(k,v)

    ss1 = mmps.Stream()
    ss1.import_file("H3O.dump", "dumppos")
    ss1.sdat.particles["pos"] -= ss1.sdat.particles["pos"][0]
    ss_tmp = copy.deepcopy(ss1)
    flag_trimming_H2O = [True for _ in range(data.total_particle)]
    H3O_data = mmps.Stream().sdat
    counter = num_of_deleted_hydrogen

    picked_idx = set()
    while(counter > 0):
        rand_idx = randrange(0, len(mol_set))
        rand_mol = mol_set[rand_idx]
        # print(rand_idx,rand_mol,dict_mol_id[rand_mol],dict_mol_type[rand_mol])
        if rand_idx in picked_idx:
            continue
        else:
            picked_idx.add(rand_idx)
            if dict_mol_type[rand_mol] != [3, 2, 2]:
                continue
            else:
                pos_original_H2O = dict_mol_pos[rand_mol]
                id_to_deleteH2O = dict_mol_id[rand_mol]
                for idx in id_to_deleteH2O:
                    flag_trimming_H2O[idx-1] = False

                axis = np.array(
                    [randrange(0, 2), randrange(0, 2), 1], dtype=float)
                angle = float(randrange(0, 180))
                rotted_pos = Rotate_particles(ss1.sdat.particles, axis, angle)
                ss_tmp.sdat.particles["pos"] = rotted_pos + pos_original_H2O[0]
                H3O_data.concate_particles(ss_tmp.sdat.particles)
                counter -= 1
                # if counter % 100 == 0:
                # print("leftParticleNum_toChangeH2OtoH3O:{:3d}".format(counter))
            counter -= 1

    return H3O_data, flag_trimming_H2O


if __name__ == '__main__':
    args = sys.argv
    inputfile = args[1]
    ss.import_file(inputfile, 'dumppos')
    inputfile = args[2]
    ss.import_file(inputfile, 'dumpbond')
    ss.sdat.set_elem_to_type(mmps.get_elem_to_type("para.rd"))
    create_molnum(ss.sdat)
    ss.sdat.sort_particles_by_id()
    print("# Success : file import has completed !")  # ファイル読み込み完了
    print("calculating now ...")

    num_of_deleted_hydrogen = delete_H_of_Nafion(ss.sdat, percentage=70)
    H3O_data, flag_trimming_H2O = change_H2O_to_H3O(
        ss.sdat, num_of_deleted_hydrogen)

    set_grid(ss.sdat, grid_len)
    latoms = get_local_atoms(ss.sdat)
    latoms = fill_o2_outside_limit_from_center(ss.sdat, latoms)
    latoms = fill_inside_center(ss.sdat, latoms, CB_radii)

    # #debug_show(ss.sdat,latoms)
    O2_data = put_O2(ss.sdat, latoms)
    # only_nafion_data = separate_NafionData(ss.sdat)
    ss.sdat.trimming_particles(flag_trimming_H2O)
    ss.sdat.concate_particles(O2_data.particles)
    ss.sdat.concate_particles(H3O_data.particles)
    
    ss.sdat.wrap_particles()

    print("number of deleted hydrogen(<- SO3H of Nafion):{:3d}".format(num_of_deleted_hydrogen))
    ss.output_file("withO2andH3O_input.rd", "input")
    ss.output_file("tmp.pos", "dumppos", ["id", "type", "pos", "velo", "force", "q", "mask", "mol"])
