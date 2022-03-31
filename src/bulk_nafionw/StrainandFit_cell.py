import copy
import sys
from MolCop import mmpystream as mmps

ss = mmps.Stream()
ifn=sys.argv[1] 
ss.import_file(ifn,'dumppos')
#print(ss.sdat.particles)

origin_cell=[u-d for u,d in zip(ss.sdat.cell,ss.sdat.dcell)]
#if you don't need to change cell, please input 0
aim_celllength =[205.127,0,0]
#change cell_length
strain = []
for idx,c in enumerate(aim_celllength):
    if c != 0:
        strain.append(c/origin_cell[idx]-1)
    else:
        strain.append(0)

ss.sdat.strain_cell(strain)
print("strain:{}".format(strain))

ss.sdat.particles['pos']+=strain

ss.output_file('fitted_'+ifn, 'dumppos', ["id", "type", "pos"])
