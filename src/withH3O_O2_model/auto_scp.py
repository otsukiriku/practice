from MolCop import mmpystream as mmps
import sys
import os
from multiprocessing import Process
import subprocess
import pexpect

#***********************************
tar_name = "dump.bond."
start_ = 0
end_ = 1000000
step_ = 100
USER = "rotsuki"
SERVER = "super.sc.imr.tohoku.ac.jp"
dir_ = "rotsuki@super.sc.imr.tohoku.ac.jp:~/rotsuki/netwark_carbon/withO2_H3O_fromver1_noweidth/"
PASSWORD = "Otumami0522!"
#***********************************

def run_scp_parallel(procs,start,end,step):
    total_procs = (end - start)//step
    job_per_1procs = total_procs//procs +1
    cnt = 0
    steps =[]
    for i in range(start,end+1,step):
        if cnt ==0:
            tmp = [i]
        cnt+=1
        if cnt == job_per_1procs:
            tmp.append(i)
            cnt=0
            steps.append(tmp)
        if i >= end:
            tmp.append(i)
            cnt=0
            steps.append(tmp)
    
    process = [0]*len(steps)
    for idx,stp in enumerate(steps):
        #print(stp[0],stp[1],step)
        process[idx] = Process(target=auto_scp, args = (stp[0],stp[1]))

    for p in process:
        p.start()
        #print(p)
    for p in process:
        p.join()
    #processes=[]
    return 0


def auto_scp(input_start,input_end):
    for i in range(input_start,input_end+1,step_):
        tar_file_name =  tar_name + str(i)
        if os.path.isfile(tar_file_name):
            #print(tar_file_name + "exist")
            continue
        target = dir_ + tar_file_name
        scp = pexpect.spawn('scp %s .' % (target))
        scp.expect('Enter*')
        scp.sendline(PASSWORD)
        scp.expect('.*ssword:')
        scp.sendline(PASSWORD)
        scp.interact()

if __name__ == "__main__":
    #procs_ = int(sys.argv[1])
    #run_scp_parallel(procs_,start_,end_,step_)
    #並列化意味ない
    auto_scp(start_,end_)
    #PASSWORD = str(input("PASSWORD?:").split())
    #print(PASSWORD)

