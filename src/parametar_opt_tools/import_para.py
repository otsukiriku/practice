import sys

def import_reaxpara(ifn):
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


if __name__ == "__main__":
    reaxpara_dict = import_reaxpara(sys.argv[1])
    print(reaxpara_dict)