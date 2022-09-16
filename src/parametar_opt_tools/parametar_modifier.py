import sys

def import_para(ifn):
    cnt = 0
    # #reaxpara dict
    reaxpara = dict()
    
    with open(ifn, 'r') as ifp:
        lines = [ l for l in ifp.readlines()]
        #header of reaxpara file
        reaxpara["header"] = lines[cnt]
        cnt+=1
        
        #start of global param
        split_line = lines[cnt].split()
        cnt+=1
        global_len = int(split_line[0])
        reaxpara["general_val"] = []
        reaxpara["general_category"] = []
        for i in range(global_len):
            split_line = lines[cnt].split()
            reaxpara["general_val"].append(float(split_line[0]))
            reaxpara["general_category"].append(split_line[1:])
            cnt += 1
        #end of global param

        #start of atom
        reaxpara["atom"] = []
        split_line = lines[cnt].split()
        atom_num = int(split_line[0])
        cnt+=4
        for a_n in range(atom_num):
            for i in range(4):
                split_line = lines[cnt].split()
                if i == 0:
                    atom_type = split_line[0]
                    tmp = [atom_type]
                    tmp.extend([ float(i) for i in split_line[1:]])
                else:
                    tmp.extend([ float(i) for i in split_line])
                cnt += 1 
            reaxpara["atom"].append(tmp)
        #end of atom

        #start of bond
        reaxpara["bond"] = []
        split_line = lines[cnt].split()
        bond_num = int(split_line[0])
        cnt += 2
        for i in range(2*bond_num):
            if i%2==0:
                split_line = lines[cnt].split()
                tmp = [ int(i) for i in split_line[:2]]
                tmp.extend([ float(i) for i in split_line[2:]])
                cnt += 1
            else:
                split_line = lines[cnt].split()
                tmp.extend([ float(i) for i in split_line])
                reaxpara["bond"].append(tmp)
                cnt += 1
        #end of bond

        #start of offdiagonal
        reaxpara["offdiagonal"] = []
        split_line = lines[cnt].split()
        offdiagonal_num = int(split_line[0])
        cnt += 1
        for i in range(offdiagonal_num):
            split_line = lines[cnt].split()
            tmp = [ int(i) for i in split_line[:2]]
            tmp.extend([ float(i) for i in split_line[2:]])
            reaxpara["offdiagonal"].append(tmp)
            cnt += 1
        #end of offdiagonal

        #start of angle
        reaxpara["angle"] = []
        split_line = lines[cnt].split()
        angle_num = int(split_line[0])
        cnt += 1
        for i in range(angle_num):
            split_line = lines[cnt].split()
            tmp = [ int(i) for i in split_line[:3]]
            tmp.extend([ float(i) for i in split_line[3:]])
            reaxpara["angle"].append(tmp)
            cnt += 1
        #end of angle

        
        #start of torsion
        reaxpara["torsion"] = []
        split_line = lines[cnt].split()
        torsion_num = int(split_line[0])
        cnt += 1
        for i in range(torsion_num):
            split_line = lines[cnt].split()
            tmp = [ int(i) for i in split_line[:4]]
            tmp.extend([ float(i) for i in split_line[4:]])
            reaxpara["torsion"].append(tmp)
            cnt += 1
        #end of torsion
        
        #start of hydrogenbond
        reaxpara["hydrogenbond"] = []
        split_line = lines[cnt].split()
        hydrogenbond_num = int(split_line[0])
        cnt += 1
        for i in range(hydrogenbond_num):
            split_line = lines[cnt].split()
            tmp = [ int(i) for i in split_line[:3]]
            tmp.extend([ float(i) for i in split_line[3:]])
            reaxpara["hydrogenbond"].append(tmp)
            cnt += 1
        #end of hydrogenbond

    return reaxpara



def output_para(reaxpara : dict, ofn = "newpara.rd"):
    header_line = reaxpara["header"]
    
    #start of general line
    general_line = "{:>3d}        ! Number of general parameters\n".format(len(reaxpara["general_val"]))
    for val,cat in zip(reaxpara["general_val"],reaxpara["general_category"]):
        general_line += "{:>9.4f} {}\n".format(val," ".join(cat))
    #end of general line

    #start of atom line
    atom_line = "{:>3d}    ! Nr of atoms; cov.r; valency;a.m;Rvdw;Evdw;gammaEEM;cov.r2;i#\n".format(len(reaxpara["atom"]))
    atom_line += "            alfa;gammavdW;valency;Eunder;Eover;chiEEM;etaEEM;n.u.\n"
    atom_line += "            cov r3;Elp;Heat inc.;n.u.;n.u.;n.u.;n.u.\n"
    atom_line += "            ov/un;val1;n.u.;val3,vval4\n"

    for each_atom_val_line in reaxpara["atom"]:
        for idx,val in enumerate(each_atom_val_line):
            if idx == 0:
                tmp = "{:>2s} ".format(val)
            elif idx == 32:
                tmp += "{:>9.4f}\n".format(val)
            elif idx % 8 == 0:
                tmp += "{:>9.4f}\n   ".format(val)
            else:
                tmp += "{:>9.4f}".format(val)
        atom_line += tmp
    #end of atom line

    #start of bond line
    bond_line = "{:>3d}      ! Nr of bonds; Edis1;LPpen;n.u.;pbe1;pbo5;13corr;pbo6\n".format(len(reaxpara["bond"]))
    bond_line += "                         pbe2;pbo3;pbo4;n.u.;pbo1;pbo2;ovcorr\n"

    for each_bond_val_line in reaxpara["bond"]:
        for idx,val in enumerate(each_bond_val_line):
            if idx == 0:
                tmp = "{:>3d}".format(val)
            elif idx == 1:
                tmp += "{:>3d}".format(val)
            elif idx == 9:
                tmp += "{:>9.4f}\n      ".format(val)
            elif idx == 17:
                tmp += "{:>9.4f}\n".format(val)
            else:
                tmp += "{:>9.4f}".format(val)
        bond_line += tmp
    #end of bond line

    #start of offdiagonal line
    offdiagonal_line = "{:>3d}    ! Nr of off-diagonal terms; Ediss;Ro;gamma;rsigma;rpi;rpi2\n".format(len(reaxpara["offdiagonal"]))
    for each_offdiagonal_val_line in reaxpara["offdiagonal"]:
        for idx,val in enumerate(each_offdiagonal_val_line):
            if idx == 0:
                tmp = "{:>3d}".format(val)
            elif idx == 1:
                tmp += "{:>3d}".format(val)
            elif idx == 7:
                tmp += "{:>9.4f}\n".format(val)
            else:
                tmp += "{:>9.4f}".format(val)
        offdiagonal_line += tmp
    #end of offdiagonal line

    #start of angle line
    angle_line = "{:>3d}    ! Nr of angles;at1;at2;at3;Thetao,o;ka;kb;pv1;pv2\n".format(len(reaxpara["angle"]))
    for each_angle_val_line in reaxpara["angle"]:
        for idx,val in enumerate(each_angle_val_line):
            if idx == 0:
                tmp = "{:>3d}".format(val)
            elif idx == 1 or idx == 2:
                tmp += "{:>3d}".format(val)
            elif idx == 9:
                tmp += "{:>9.4f}\n".format(val)
            else:
                tmp += "{:>9.4f}".format(val)
        angle_line += tmp
    #end of angle line

    #start of torsion line
    torsion_line = "{:>3d}    ! Nr of torsions;at1;at2;at3;at4;;V1;V2;V3;V2(BO);vconj;n.u;n\n".format(len(reaxpara["torsion"]))
    for each_torsion_val_line in reaxpara["torsion"]:
        for idx,val in enumerate(each_torsion_val_line):
            if idx == 0:
                tmp = "{:>3d}".format(val)
            elif idx == 1 or idx == 2 or idx == 3:
                tmp += "{:>3d}".format(val)
            elif idx == 10:
                tmp += "{:>9.4f}\n".format(val)
            else:
                tmp += "{:>9.4f}".format(val)
        torsion_line += tmp
    #end of torsion line

    #start of hydrogenbond line
    hydrogenbond_line = "{:>3d}    ! Nr of hydrogen bonds;at1;at2;at3;Rhb;Dehb;vhb1\n".format(len(reaxpara["hydrogenbond"]))
    for each_hydrogenbond_val_line in reaxpara["hydrogenbond"]:
        for idx,val in enumerate(each_hydrogenbond_val_line):
            if idx == 0:
                tmp = "{:>3d}".format(val)
            elif idx == 1 or idx == 2:
                tmp += "{:>3d}".format(val)
            elif idx == 6:
                tmp += "{:>9.4f}\n".format(val)
            else:
                tmp += "{:>9.4f}".format(val)
        hydrogenbond_line += tmp
    #end of hydrogenbond line

    with open(ofn, 'w') as ofp:
        ofp.write(header_line+general_line+atom_line+bond_line+offdiagonal_line+angle_line+torsion_line+hydrogenbond_line)



def modify_reaxpara_relatively(reaxpara : dict,key, target, idx, differ):
        """
        reaxffのパラメーターを調整したいときに使う
        :param key:reaxparaの操作したい項目
        :param target:操作したい原子種のリスト
        :param idx:変更したい項目のインデックス(※1スタートのインデックス)
        :param differ:変更したい差分の値
        """
        if idx <=0:
            print("Error : index should start with 1")
            return
        idx -= 1
           
        if key == "general_val":
            reaxpara["general_val"][idx] += differ
        elif key == "atom":
            for val in reaxpara["atom"]:
                if target == val[0]:
                    val[idx] += differ
        elif key == "bond":
            for val in reaxpara["bond"]:
                if target == val[:2]:
                    val[idx] += differ
        elif key == "offdiagonal":
            for val in reaxpara["offdiagonal"]:
                if target == val[:2]:
                    val[idx] += differ
        elif key == "angle":
            for val in reaxpara["angle"]:
                if target == val[:3]:
                    val[idx] += differ
        elif key == "torsion":
            for val in reaxpara["torsion"]:
                if target == val[:4]:
                    val[idx] += differ
        elif key == "hydrogenbond":
            for val in reaxpara["hydrogenbond"]:
                if target == val[:3]:
                    val[idx] += differ
        else:
            print("Error : key was not founded in parameter")
            print("Available key :  \"general_val\", \"atom\", \"bond\", \"offdiagonal\",\"angle\",\"torsion\",\"hydrogenbond\"")

