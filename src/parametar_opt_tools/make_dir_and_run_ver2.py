import sys
import os
import shutil
import subprocess
import math
from E_T.func import parametar_modifier as para_modi
"""
template ディレクトリが同一ディレクトリにないと動かない．
"""

#********************************************
PATH_TO_TEMPLATE_DIR = "./template"
TAR_CATEGORY = "angle"
TAR_NAME = "pval1"
TAR_LIST_OF_ATOM = [2,2,3]
TAR_IDX_1START_IN_ARRAY = 5
#********************************************


def run_auto_run():
    #result = subprocess.run(["python3","../template/only_analysis_autorun.py"])
    result = subprocess.run(["python3","auto_run.py"])
    return result

def main(min_range_l,max_range_l,step_l):
    #org_param_val = get_org_param_val()
    os.chdir(PATH_TO_TEMPLATE_DIR)
    paradict = para_modi.import_para("para0726.rd")
    org_param_val = para_modi.get_reaxpara(paradict,TAR_CATEGORY,TAR_LIST_OF_ATOM,TAR_IDX_1START_IN_ARRAY)
    os.chdir("../")
    step_num = int(math.ceil((max_range_l - min_range_l)/step_l)) + 1
    range_lst = [ min_range_l+round(i*step_l,4) for i in range(step_num)]
    for i in range_lst:
        if i == 0.0:
            continue
        new_dir_name = TAR_NAME +"_"+ str(round(org_param_val+i,4))
        print(new_dir_name)
        if not os.path.exists(new_dir_name):
            shutil.copytree(PATH_TO_TEMPLATE_DIR,"./"+new_dir_name)
            os.chdir(new_dir_name)
            paradict = para_modi.import_para("para0726.rd")
            new_paradict = para_modi.modify_reaxpara_relatively(paradict,TAR_CATEGORY,TAR_LIST_OF_ATOM,TAR_IDX_1START_IN_ARRAY,i)
            para_modi.output_para(paradict,"para0726.rd")
            os.chdir("../")
        os.chdir(new_dir_name)
        run_auto_run()
        os.chdir("../")
if __name__ == "__main__":
    min_range = sys.argv[1]
    if min_range.isdecimal():
        min_range = int(min_range)
    else:
        min_range = float(min_range)
    max_range = sys.argv[2]
    if max_range.isdecimal():
        max_range = int(max_range)
    else:
        max_range = float(max_range)
    step = sys.argv[3]
    if step.isdecimal():
        step = int(step)
    else:
        step = float(step)
    main(min_range,max_range,step)
