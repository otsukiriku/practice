from MolCop import mmpystream as mmps
from MolCop.analysis import topology
from E_T.func import Pt_location
from E_T.func import read_center as rc
from E_T.func import get_center as gc
import matplotlib.pyplot as plt
import numpy as np
#import pandas as pd
from sklearn.linear_model import LinearRegression
import sys
import copy
import math
sys.path.append("/nfshome12/rotsuki/molcop/src/")
#print(sys.path)
import modeling.unwrap_particles as u_p
import evaluate_structure.get_center as g_c
sys.path.append("/nfshome12/rotsuki/practice/analysis_src/")
import  cover_analysis_bymask as c_a

#########################################
Pt_cluster_num=309
Pt_num=56
Pt_mask_start = 7
#cut off distance from Pt surface
cutoff = 9
target_type = 6
dot_color=["blue","red", "green", "orange","b"]
line_color=["c","pink","olive","y","b"]
#########################################

def plot_scatter(xdata:list, ydata:list, _figname = "test", xlabel='x', ylabel='y'):
    """datax,y should be 2dimension list"""
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    #ax.set_ylim(0, 10)
    ax.set_xlabel("Prove Radius(Aungstrom)", size = 16,)
    ax.set_ylabel("Coverage(%)", size = 16,)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    ax.set_xlim(5,45)
    ax.set_ylim(0,80)
    #plt.figure(figsize=(4,3))
    print("xdata dimension")
    print(len(xdata))
    if len(xdata) == 1:
        ax.scatter(xdata,ydata)
        test = [[i] for i in xdata]
        X = np.array(test)
        Y = np.array(ydata)
        Ir = LinearRegression()
        Ir.fit(X,Y)
        plt.plot(X,Ir.predict(X), color="red")
        ax.text(0.95,0.05,"coefficient="+str(Ir.coef_), transform=ax.transAxes, horizontalalignment="right")
        fig.savefig(_figname+".png",dpi=250)
        plt.show()
    else:
        for idx, (x,y) in enumerate(zip(xdata,ydata)):
            ax.scatter(x,y,c=dot_color[idx])
            test = [[i] for i in x]
            X = np.array(test)
            Y = np.array(y)
            Ir = LinearRegression()
            Ir.fit(X,Y)
            plt.plot(X,Ir.predict(X), color=line_color[idx])
            ax.text(0.95,0.05+0.05*idx, dot_color[idx]+":coefficient="+str(Ir.coef_), transform=ax.transAxes, horizontalalignment="right")
        fig.savefig(_figname+".png",dpi=250)

def main(ss:mmps.Stream()):
    CB_center = rc.r_c("./center.txt")

    mask, prove = Pt_location.Pt_location(ss.sdat, nonzero_option=False, prove_range = 40)

    ss.sdat.create_connect_list(0.3)
    connect_list = ss.sdat.connect_list

    #flagconnect reindexの両方をTrueにすることで
    #connect_listを残してトリミングすることもできる．
    Pt = copy.deepcopy(ss.sdat)
    flag = ss.sdat.particles['type']==4
    Pt.flagconnect = True
    Pt.trimming_particles(flag, reindex=True)

    #配位数でPtの表面を管理することにした．
    coordination_num = c_a.get_coordination(Pt)
    distflag = c_a.distance_flag(Pt, CB_center)
    flag = (coordination_num < 10) & distflag
    Pt.trimming_particles(flag, reindex=True)
    topology.unwrap_molecule(Pt)
    m_list = topology.create_molecule(Pt)
    Pt.wrap_particles()
    #Ptはこれで表面だけ残った．
    #trim only target particle
    target_flag = ss.sdat.particles['type'] == target_type
    ss.sdat.trimming_particles(target_flag)

    target = ss.sdat

    #print(cover)
    Pt_mask = [i for i in range(Pt_mask_start, Pt_mask_start+Pt_num)]
    count, cover = c_a.calc_coverage(target, Pt, Pt_mask)
    #print("count")
    #print(count)
    #print("cover")
    #print(cover)
    print('Pt_mask Prove Number_of_target_atom coverage(%)')

    coverage_arr =[]
    for idx, c in enumerate(count):
        total_Pt = len(Pt.particles['id'][Pt.particles['mask']==idx+Pt_mask_start])
        coverage = cover[idx]/total_Pt*100
        #covered_Pt = np.sum(Pt.particles['flag'][Pt.particles['mol']==idx+1])
        print("{} {:.3} {} {:.3}%".format(idx+Pt_mask_start, prove[idx], c, coverage))
        coverage_arr.append(coverage)
    non0coverage = []
    non0prove = []
    for cov,prov in zip(coverage_arr,prove):
        if prov != 0.0:
            non0coverage.append(cov)
            non0prove.append(prov)
    #print(non0prove)
    return non0coverage,non0prove

    """
    # for debug trimmed_Pt
    ss1 = mmps.Stream()
    ss1.sdat = Pt
    ss1.output_file("show_dump", 'dumppos', ["id","type","pos",'mol','mask'] )
    ss.output_file("trimmed_input", 'input')
    """

if __name__ == "__main__":
    #for pandas
    #data_num=int(len(sys.argv)/2)
    #colum_rabel=[ "data"+str(i) for i in range(0,data_num)]

    ss = mmps.Stream()
    ss.import_file(sys.argv[1], 'dumppos')
    ss.import_file(sys.argv[2], 'dumpbond')
    ss_coverage, ss_prove=main(ss)
    #ss_coverage = np.array([ss_coverage])
    #ss_prove = np.array([ss_prove])
    print("ss_coverage")
    print(ss_coverage)
    ss1 = mmps.Stream()
    ss1.import_file(sys.argv[3], 'dumppos')
    ss1.import_file(sys.argv[4], 'dumpbond')
    ss1_coverage, ss1_prove = main(ss1)
    print("ss1_coverage")
    print(ss1_coverage)
    ss_coverage=[ss_coverage]
    ss_coverage.append(ss1_coverage)
    ss_prove=[ss_prove]
    ss_prove.append(ss1_prove)
    #ss_coverage=np.append(ss_coverage,np.array([ss1_coverage]),axis=0)
    #ss_prove=np.append(ss_prove,np.array([ss1_prove]),axis=0)
    print("ss_coverage")
    print(ss_coverage)
    #print("ss_prove")
    #print(ss_prove)
    plot_scatter(ss_prove, ss_coverage)

