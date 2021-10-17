#!/usr/bin/env python3

#!/usr/bin/env python3

"""
How to use

execution
→file name(split names with space)
e.g.)



"""



import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
filename = input().split()
print(filename)

dfname = []

for fp in filename:

    file = open(fp, mode='r')
    a=file.readline().rstrip('\n')

    if a == '# Radial distribution function (200 data points):' :
        title = file.readline().rstrip('\n')
    else:
        title = a

    title=title.replace('# "Pair separation distance" ','')
    title=title.split(' ')
    title.append('none')
    title.insert(0,"x")
    #print(title)

    dfn = pd.read_csv(fp, sep = " ",skiprows=2,names=title)
    dfname.append(dfn)
    file.close()
    #print(dfn.head())
    #print(dfn.keys)


colors=['b','r','b','g']
#withOH
#O of H2O
fig = plt.figure(figsize=(4, 3))
ax = fig.add_subplot(111)
ax.set_ylim(0, 40)
ax.set_xlabel("r(Å)", size = 16,)
ax.set_ylabel('g(r)', size = 16,)
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
for dfn,cl in zip(dfname,colors):
    ax.plot(dfn['x'],dfn['Pt-OofH2O'],color = cl)
fig.savefig("Pt-OofH2O"+"multiple"+".png",dpi=250)


#C of Nafion
fig = plt.figure(figsize=(4, 3))
ax = fig.add_subplot(111)
ax.set_ylim(0, 2)
ax.set_yticks([0,0.5,1.0,1.5,2.0])
ax.set_xlabel("r(Å)", size = 16,)
ax.set_ylabel('g(r)', size = 16,)
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
for dfn,cl in zip(dfname,colors):
    ax.plot(dfn['x'],dfn['CofNafion-Pt'],color = cl)
fig.savefig("CofNafion-Pt"+"multiple"+".png",dpi=250)


#S of Nafion
fig = plt.figure(figsize=(4, 3))
ax = fig.add_subplot(111)
ax.set_ylim(0, 5)
ax.set_yticks([0,0.5,1.0,1.5,2.0])
ax.set_xlabel("r(Å)", size = 16,)
ax.set_ylabel('g(r)', size = 16,)
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
for dfn,cl in zip(dfname,colors):
    ax.plot(dfn['x'],dfn['Pt-SofNafion'],color = cl)
fig.savefig("SofNafion-Pt"+"multiple"+".png",dpi=250)
