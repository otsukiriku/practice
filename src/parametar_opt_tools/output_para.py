import sys

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


if __name__ == "__main__":
    reaxpara_dict = output_para(sys.argv[1])
    print(reaxpara_dict)