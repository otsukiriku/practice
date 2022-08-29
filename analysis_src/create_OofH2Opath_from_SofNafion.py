#!/usr/bin/env python3
from MolCop import mmpystream as mmps
from MolCop.analysis import topology
from MolCop.analysis import neighbor
import collections
import numpy as np
import sys
import copy
import math

#*********************************************************
#cut off distance from Pt surface
cluster_type = 13
cutoff = 3.5
start_search_type = 5
# 5:SofNafion 6:FofNadion
#*********************************************************

#ss = mmps.Stream()
#ss_local = copy.deepcopy(ss.sdat)
def find_cluster_from(ss_local:mmps.Stream().sdat, start_type=5, find_type=13,CL=3.5, ):
    #flag = (ss_local.particles["type"] == cluster_type) | (ss_local.particles["type"] == start_search_type)
    start_flag = ss_local.particles["type"] == start_type
    find_flag = ss_local.particles["type"] == find_type
    n_list = neighbor.make_neighbor(ss_local,cutoff)

    visited = [0 for _ in range(ss_local.total_particle)]
    start_idx = ss_local.particles["id"][start_flag]
    for s in start_idx:
        visited[s-1] = 1
        que = collections.deque()
        que.append(n_list[s-1])
        while que:
            tmp = que.popleft()
            #print(tmp)
            for t in tmp:
                if visited[t]==0 and find_flag[t]:
                    que.append(n_list[t])
                    visited[t]=1
                else:
                    continue

    if "cluster_flag" not in ss_local.particles:
        ss_local.add_particles_property("cluster_flag", _dtype=int, dim=1)
    ss_local.particles["cluster_flag"]=np.array(visited)
    return ss_local

if __name__ == "__main__":
    ss_out = mmps.Stream()
    ss_out.import_file(sys.argv[1], "dumppos")
    ss_out.sdat = find_cluster_from(ss_out.sdat,start_search_type,cluster_type,cutoff)
    ss_out.output_file("withclusterinfo_"+sys.argv[1],"dumppos",["id","type","pos","cluster_flag"])
    
