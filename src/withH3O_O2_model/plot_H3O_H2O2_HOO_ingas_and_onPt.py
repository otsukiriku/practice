
from MolCop import mmpystream as mmps
#from MolCop.analysis import topology
#from E_T.func import figure
#from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import numpy as np
import sys
import math
import os

dot_color=["black","salmon","red","lightblue","blue", "lightgreen","green"]

def plot_scatter(xdata:list, ydata:list, _figname = "test", xlabel='x', ylabel='y',_dot_label=None):
    ydata = np.array(ydata)
    if ydata.ndim == 1:
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)
        #ax.set_ylim(0, 10)
        #ax.set_ylim(0, 3500)
        ax.set_xlabel(xlabel, size = 16,)
        ax.set_ylabel(ylabel, size = 16,)
        ax.tick_params(axis='x', labelsize=16)
        ax.tick_params(axis='y', labelsize=16)
        #plt.figure(figsize=(4,3))
        ax.scatter(xdata,ydata)
        fig.savefig(_figname+".png",dpi=250)
    else:
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)
        #ax.set_ylim(0, 10)
        #ax.set_ylim(0, 3500)
        ax.set_xlabel(xlabel, size = 16,)
        ax.set_ylabel(ylabel, size = 16,)
        ax.tick_params(axis='x', labelsize=16)
        ax.tick_params(axis='y', labelsize=16)
        for idx,y in enumerate(ydata):
            ax.scatter(xdata,y,color=dot_color[idx])
            ax.text(0.50,0.95-0.06*idx, dot_color[idx]+":"+_dot_label[idx], transform=ax.transAxes, horizontalalignment="left", fontsize=11)
        fig.savefig(_figname+".png",dpi=250)

def main(start_step,end_step, _dirname=None):
    ss = mmps.Stream()
    ss.import_file('../calc/config.rd', 'config')

    total_step = end_step - start_step

    timestep = []
    isH3O = []
    isH2O2_in_gas = []
    isH2O2_on_Pt = []
    isHOO_in_gas = []
    isHOO_on_Pt = []
    isH2O_in_gas_from_O2gas = []
    isH2O_on_Pt_from_O2gas = []
    isO2gas = []

    for step in range(start_step, end_step+1, ss.sdat.file_step):
    #for step in range(start_step, end_step+1, 1000):
        ss = mmps.Stream()
        ifn_pos = '../color/with_H3O_H2O2_HOO_search_colored.pos.' +  str(step)
        if not os.path.isfile(ifn_pos):
            print(f'{ifn_pos} is not exit')
            continue
        ss.import_file(ifn_pos, 'dumppos')

        timestep.append(step)
        isH3O.append(np.sum(ss.sdat.particles["isH3O_in_gas"])//4)
        isH2O2_in_gas.append(np.sum(ss.sdat.particles["isH2O2_in_gas"])//4)
        isH2O2_on_Pt.append(np.sum(ss.sdat.particles["isH2O2_on_Pt"])//4)
        isHOO_in_gas.append(np.sum(ss.sdat.particles["isHOO_in_gas"])//3)
        isHOO_on_Pt.append(np.sum(ss.sdat.particles["isHOO_on_Pt"])//3)
        isH2O_in_gas_from_O2gas.append(np.sum(ss.sdat.particles["isH2O_in_gas_from_O2gas"])//3)
        isH2O_on_Pt_from_O2gas.append(np.sum(ss.sdat.particles["isH2O_on_Pt_from_O2gas"])//3)
        #isO2gas.append(np.sum(ss.sdat.particles["isO2gas"])//2)

    os.chdir("../")
    if os.path.exists("result"):
        os.chdir("result")
    else:
        os.mkdir("result")
        os.chdir("result")

    if _dirname is None:
        plot_scatter(timestep,isH3O,"H3Onum_vs_timestep",xlabel="timestep",ylabel="isH3O")
        plot_scatter(timestep,isH2O2_in_gas,"H2O2_in_gas_vs_timestep",xlabel="timestep",ylabel="isH2O2_in_gas")
        plot_scatter(timestep,isH2O2_on_Pt,"H2O2_on_Pt_vs_timestep",xlabel="timestep",ylabel="isH2O2_on_Pt")
        plot_scatter(timestep,isHOO_in_gas,"HOO_in_gas_vs_timestep",xlabel="timestep",ylabel="isHOO_in_gas")
        plot_scatter(timestep,isHOO_on_Pt,"HOO_on_Pt_vs_timestep",xlabel="timestep",ylabel="isHOO_on_Pt")
        plot_scatter(timestep,isH2O_in_gas_from_O2gas,"H2O_in_gas_from_O2gas_vs_timestep",xlabel="timestep",ylabel="isH2O_in_gas_from_O2gas")
        plot_scatter(timestep,isH2O_on_Pt_from_O2gas,"H2O_on_Pt_from_O2gas_vs_timestep",xlabel="timestep",ylabel="isH2O_on_Pt_from_O2gas")
        plot_scatter(timestep,[isH3O,isH2O2_in_gas,isH2O2_on_Pt,isHOO_in_gas,isHOO_on_Pt,isH2O_in_gas_from_O2gas,isH2O_on_Pt_from_O2gas],"H3O,H2O2,HOOnum_vs_timestep",xlabel="timestep",ylabel="isH3O,isH2O2,isHOO,isH2O",_dot_label=["isH3O","isH2O2_in_gas","isH2O2_on_Pt","isHOO_in_gas","isHOO_on_Pt","isH2O_in_gas_from_O2gas","isH2O_on_Pt_from_O2gas"])
    else:
        plot_scatter(timestep,isH3O,_dirname+"H3Onum_vs_timestep",xlabel="timestep",ylabel="isH3O")
        plot_scatter(timestep,isH2O2_in_gas,_dirname+"H2O2_in_gas_vs_timestep",xlabel="timestep",ylabel="isH2O2_in_gas")
        plot_scatter(timestep,isH2O2_on_Pt,_dirname+"H2O2_on_Pt_vs_timestep",xlabel="timestep",ylabel="isH2O2_on_Pt")
        plot_scatter(timestep,isHOO_in_gas,_dirname+"HOO_in_gas_vs_timestep",xlabel="timestep",ylabel="isHOO_in_gas")
        plot_scatter(timestep,isHOO_on_Pt,_dirname+"HOO_on_Pt_vs_timestep",xlabel="timestep",ylabel="isHOO_on_Pt")
        plot_scatter(timestep,isH2O_in_gas_from_O2gas,_dirname+"H2O_in_gas_from_O2gas_vs_timestep",xlabel="timestep",ylabel="isH2O_in_gas_from_O2gas")
        plot_scatter(timestep,isH2O_on_Pt_from_O2gas,_dirname+"H2O_on_Pt_from_O2gas_vs_timestep",xlabel="timestep",ylabel="isH2O_on_Pt_from_O2gas")
        #plot_scatter(timestep,isO2gas,_dirname+"O2gas_vs_timestep",xlabel="timestep",ylabel="isO2gas_from_O2gas")
        plot_scatter(timestep,[isH3O,isH2O2_in_gas,isH2O2_on_Pt,isHOO_in_gas,isHOO_on_Pt,isH2O_in_gas_from_O2gas,isH2O_on_Pt_from_O2gas],_dirname+"H3O,H2O2,HOOnum_vs_timestep",xlabel="timestep",ylabel="isH3O,isH2O2,isHOO,isH2O",_dot_label=["isH3O","isH2O2_in_gas","isH2O2_on_Pt","isHOO_in_gas","isHOO_on_Pt","isH2O_in_gas_from_O2gas","isH2O_on_Pt_from_O2gas"])

if __name__ == "__main__":
    s_step = int(sys.argv[1])
    e_step = int(sys.argv[2])
    if len(sys.argv)>3:
        name = sys.argv[3]
        main(s_step,e_step,name)
    else:
        main(s_step,e_step)
