from MolCop import mmpystream as mmps
import sys
import numpy as np
import itertools


def create_shakepara(atoms:mmps.Stream().sdat, shake_list:list, shake_type:list, mask_start):
    para_list = []
    shake_type_set = list(set(shake_type))
    #print(shake_type_set)
    shake_mask = [i+mask_start for i in range(len(shake_type))]

    atoms.para_dict ={
        "OH":[1, 9.7],
        "H2O":[4, 0.970, 0.970, 1.525],
        "CH":[1, 9.7],
        "CH2":[2, 1.090, 1.090],
        "CH3":[3, 1.090, 1.090, 1.090]
    }

    for msk,s_t in zip(shake_mask,shake_type_set):
        atoms.para_dict[s_t].insert(0,msk)
        para_list.append(atoms.para_dict[s_t])
    atoms.set_shakepara(para_list)
    print(atoms.para_dict)

def create_shakelist(atoms:mmps.Stream().sdat, shake_list:list, shake_type:list):
    shake_type_num = [ atoms.para_dict[t][1] for t in shake_type ]
    atoms.set_shakelist(shake_type_num,shake_list)
    if not 'mask' in atoms.particles:
        atoms.add_particles_property('mask', _dtype = int, dim=1)
    for s_t,s_l in zip(shake_type,shake_list):
        s_l = [i-1 for i in s_l]
        s_l.pop(0)
        #print(s_t)
        #print(atoms.para_dict[s_t])
        atoms.particles['mask'][s_l] = atoms.para_dict[s_t][0]
        #atoms.particles['mask'] = 

def search_shakelist(atoms:mmps.Stream().sdat):
    bonds = [ [] for _ in range(len(atoms.particles['id']))]
    shake_type = [ False for _ in range(len(atoms.particles['id']))]
    #print(bonds)

    #atoms.add_particles_property('Hconnect', _dtype = int, dim=1)

    H_flag = atoms.particles['type'] == 2
    connect_with_H = []
    #Hと結合している原子のconnectlistをconnect_with_Hに格納する
    for c_list,fl in zip(atoms.connect_list,H_flag):
        if fl:
            connect_with_H.append(c_list)
    H_id = atoms.particles['id'][H_flag]
    #print(connect_with_H)
    #くっついている水素のインデックスをbondsに格納する．
    for connect_1d, H in zip(connect_with_H, H_id):
        for c in connect_1d:
            bonds[c].append(H)
    #print(bonds)

    allbond = list(itertools.chain.from_iterable(bonds))
    allbond = list(set(allbond))
    shake_type = []
    for idx,b in enumerate(bonds):
        if b:
            Hbond_num = len(b)
            root_type = atoms.particles['type'][idx]
            #Hconnect with O
            key = "NoMatch"
            if root_type == 3:
                if Hbond_num == 1:
                    key = "OH"
                if Hbond_num == 2:
                    key = "H2O"

            #Hconnect with C
            if root_type == 1:
                if Hbond_num == 1:
                    key = "CH"
                if Hbond_num == 2:
                    key = "CH2"
                if Hbond_num == 3:
                    key = "CH3"

            b.insert(0,idx+1)
            b.insert(0,key)
            shake_type.append(key)
    shake_list = [ b for idx,b in enumerate(bonds) if atoms.particles['type'][idx] != 2]
    shake_list = [ s for s in shake_list if s]
    shake_list = [ s for s in shake_list if s[0] != "NoMatch"]
    shake_type = [ s.pop(0) for s in shake_list]
    #print(shake_list)
    #print(shake_type)
    #print(allbond)
    return shake_list, shake_type

if __name__ == "__main__":
    ss = mmps.Stream()
    ss.import_file(sys.argv[1],'dumppos')
    ss.import_file(sys.argv[2],'dumpbond')

    c_list = ss.sdat.create_connect_list()
    #ss.sdat.add_particles_property()
    shake_list,shake_type = search_shakelist(ss.sdat)
    create_shakepara(ss.sdat,shake_list,shake_type,mask_start = 90)
    create_shakelist(ss.sdat,shake_list,shake_type)
    ss.output_file("show_dump", "dumppos",['id','type','pos','mask'])
    ss.output_file("test_input", "input")
