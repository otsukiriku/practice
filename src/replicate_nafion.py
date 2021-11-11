import copy
from MolCop import mmpystream as mmps

ss = mmps.Stream()
ss1 = mmps.Stream()
ss2 = mmps.Stream()
ss_out = mmps.Stream()

ss.import_file('nafionw.dump', 'dumppos')
#print(ss.sdat.particles)

flag = ss.sdat.particles['id']<124417
#print(flag)
H2O = ss1.sdat
H2O.particles["id"] = ss.sdat.particles['id'][flag]
H2O.particles["type"] = ss.sdat.particles['type'][flag]
H2O.particles["pos"] = ss.sdat.particles['pos'][flag]
#print(H2O.particles)
H2O.replicate_particles([1,1,3])
ss_out.sdat.concate_particles(H2O.particles)


flag = ss.sdat.particles['id']>124416
#print(flag)
Naf = ss2.sdat
Naf.particles["id"] = ss.sdat.particles['id'][flag]
Naf.particles["type"] = ss.sdat.particles['type'][flag]
Naf.particles["pos"] = ss.sdat.particles['pos'][flag]
#print(H2O.particles)
Naf.replicate_particles([1,1,3])
ss_out.sdat.concate_particles(Naf.particles)

print(ss_out)

ss_out.output_file('newnafionw.dump', 'dumppos', ["id", "type", "pos"])
