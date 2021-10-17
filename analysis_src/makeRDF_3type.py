#!/usr/bin/env python3
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
filename = input()


file = open(filename, mode='r')
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


df = pd.read_csv(filename, sep = " ",skiprows=2,names=title)
file.close()
#print(df.head())
#withOH
    
#O of H2O
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111)
x=df['x'].values
y=df['Pt-OofH2O'].values
ax.set_ylim(0, 40)
ax.set_xlabel("r(Å)", size = 16,)
ax.set_ylabel('g(r)', size = 16,)
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
#plt.figure(figsize=(4,3))
ax.plot(x,y)
fig.savefig("Pt-OofH2O"+filename+".png",dpi=250)
plt.show()

#C of Nafion#O of H2O
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111)
x=df['x'].values
y=df['CofNafion-Pt'].values
ax.set_ylim(0, 2)
ax.set_yticks([0,0.5,1.0,1.5,2.0])
ax.set_xlabel("r(Å)", size = 16,)
ax.set_ylabel('g(r)', size = 16,)
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
#plt.figure(figsize=(4,3))
ax.plot(x,y)
fig.savefig("Pt-CofNafion"+filename+".png",dpi=250)
plt.show()

#S of Nafion
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111)
x=df['x'].values
y=df['Pt-SofNafion'].values
ax.set_ylim(0, 5)
ax.set_xlabel("r(Å)", size = 16,)
ax.set_ylabel('g(r)', size = 16,)
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
#plt.figure(figsize=(4,3))
ax.plot(x,y)
fig.savefig('Pt-SofNafion'+filename+".png",dpi=250)
plt.show()
