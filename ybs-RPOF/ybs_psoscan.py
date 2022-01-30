import os
import math, re
import numpy as np
import matplotlib.pyplot as plt
#
from threading import Thread
#
from time import sleep
from ybs_ffieldpam import ffpam
from ybs_write_params import w_params
from ybs_mul_ffield import mul_ffield as m_ffx
from ybs_pso import PSO
from ybs_id_v_flag_turn import id_v_flag_turn as pigturn
def psoscan(ff,f,scale=0 , num_shunxu = 1 , num_luanxu = 20 , c1 = 1.8 , c2 = 1.5 ,w = 0.5 ,num_pso = 1 ,POP_SIZE = 8 ,pso_x_num = 6, pam_={},outfile='outdata'):

    num_pso = num_pso        #
    POP_SIZE = POP_SIZE      #
    pso_x_num = pso_x_num    #
    num_shunxu = num_shunxu  #
    num_luanxu = num_luanxu  #
    w = w                    #
    c1 = c1                  #
    c2 = c2                  #

    scale_flag = scale
    outdir = outfile
    dir_path = os.getcwd()
    ############################################################################################# def func
    def pso_scan(f, m_ffx, x_pam_start, var_name, bound,pam,f_x_flag,pso_start_x_fx_dir,num_pso=1,POP_SIZE=20):
        pattern = r',|\|/|;|\'|`|\[|\]|<|>|\?|:|"|\{|\}|\~|!|@|#|\$|%|\^|&|\(|\)|-|=|\_|\+|，|。|、|；|‘|’|【|】|·|！| |…|（|）'
        print("yes")
        f = f
        m_ffx = m_ffx
        x_pam_start = x_pam_start
        var_name = var_name
        bound = bound
        POP_SIZE = POP_SIZE
        f_x_flag = {}
        #w = 0.5
        #c1 = 1.8
        #c2 = 1.5
        v_max = bound[0][1]-bound[0][0]
        x_fx_list = {}
        for v_d in bound:
            v_dx = v_d[1]-v_d[0]
            if v_dx < v_max:
                v_max = v_dx
        if len(var_name) == 1:
            x_fx_list = pso_start_x_fx_dir
            pam_flag = var_name[0]
            for k_pam,v_fx in x_fx_list.items():
                pam_v = k_pam
                f_v = v_fx
                x_flag_value_str_list = pam_flag + "_" + str(pam_v) + "_"
                f_x_flag[x_flag_value_str_list] = f_v
        print(x_pam_start, var_name, bound)
        with open("outdata/pso_cal", 'w') as www1:
            www1.write("")
        with open("outdata/pso_result", 'w') as ww1:
            ww1.write("")
        ## def
        def func(*x_list):
            x_flag_value = {}
            x_flag_value_str_list = ""
            for i in range(len(var_name)):
                pam_flag = var_name[i]
                pam_v = float("%.4f" % x_list[i])
                x_flag_value[pam_flag] = pam_v
                x_flag_value_str = pam_flag + "_" + str(pam_v) + "_"
                x_flag_value_str_list += x_flag_value_str
            # pigturn
            x_flag_value_str_list_pre = x_flag_value_str_list
            x_flag_value_str_list = pigturn(pam,x_flag_value_str_list)[1]
            pso_x_name = m_ffx(x_flag_value, "pso")
            if x_flag_value_str_list in f_x_flag:
                f_v = -f_x_flag[x_flag_value_str_list]
            else:
                f_v = -f(pso_x_name)[0]
                f_x_flag[x_flag_value_str_list] = -f_v
                # 单
                x_fx_list[pam_v] = -f_v
                with open("outdata/pso_cal", 'a') as www:
                    www.write((x_flag_value_str_list_pre + "  ").ljust(40) + str(-f_v) + "\n")
            return f_v
        ## end def
        pso = PSO(func, bound, POP_SIZE, x_pam_start, w, c1, c2, v_max)
        pso_result = {}
        ##
        for _ in range(num_pso):
            pso_result[_] = []
            pso_result[_] = pso.pso()

        min_f_v = 99999999999999
        x_re_value = {}
        for _ in range(num_pso):
            result = pso_result[_]
            pso_re_fx_v = result[-1]
            if pso_re_fx_v < min_f_v:
                min_f_v = pso_re_fx_v
                min_x_v = result[0]
        if 1 == 2:
            with open("outdata" + os.sep + "pso_result") as r:
                lines = r.readlines()
                min_f_v = 99999999999999
                min_x_v = []
                x_re_value = {}
                for line in lines:
                    line = line.replace('[', '').replace(']', '')
                    spline = line.split()
                    if len(spline) != 0:
                        pso_re_pam_val_list = spline[:-1]
                        print(pso_re_pam_val_list)
                        pso_re_fx_v = float(spline[-1])
                        if pso_re_fx_v < min_f_v:
                            min_f_v = pso_re_fx_v
                            min_x_v = pso_re_pam_val_list
        for i in range(len(var_name)):
            pam_flag = var_name[i]
            pam_a = bound[i][0]
            pam_b = bound[i][1]
            id1 = pam_flag.split("-")[0]
            id2 = pam_flag.split("-")[1]
            id3 = pam_flag.split("-")[2]
            pam_v = float("%.4f" % float(min_x_v[i]))
            x_re_value[pam_flag] = pam_v
            with open(os.getcwd() + os.sep + outdir + os.sep + 'xkpambest', 'a') as w_pambest:
                w_pambest.write(str(id1).ljust(3) + str(id2).ljust(3) + str(id3).ljust(5) \
                                + str(pam_v).rjust(10) + str(pam_a).rjust(10) + str(pam_b).rjust(10) + '\n')
        with open(os.getcwd() + os.sep + outdir + os.sep + 'xkpambest', 'a') as w_pambest:
            w_pambest.write("#".ljust(50) + str(min_f_v) +'\n')
        pso_xf_best = m_ffx(x_re_value, "ffield.reax.start")
        return x_re_value,x_fx_list ,min_f_v

    def pam_pso_scan(pam,f_x_flag,pso_x_num = 5,shunxu=0):
        numpam = len(pam)
        pso_x_num = pso_x_num
        numpam_list = list(range(1,numpam + 1))
        if shunxu == 0:
            pam_step_list = [numpam_list[i:i+pso_x_num] for i in range(0,len(numpam_list),pso_x_num)]
        else:
            numpam_list = np.random.permutation(numpam_list)

            pam_step_list = [numpam_list[i:i + pso_x_num] for i in
                             range(0, len(numpam_list), pso_x_num)]
        for pam_step in pam_step_list:
            var_name = []
            bound = []
            x_pam_start = []
            #num_pso = 1
            #POP_SIZE = 8

            for num in pam_step:
                pam_a = pam[num]['value'][1]
                pam_b = pam[num]['value'][2]
                id1 = int(pam[num]['id'][1])
                id2 = int(pam[num]['id'][2])
                id3 = int(pam[num]['id'][3])
                pam_v = ffpam(xstart, [id1, id2, id3])
                pam_flag = str(id1) + '-' + str(id2) + '-' + str(id3)
                var_name.append(pam_flag)
                x_pam_start.append(pam_v)
                bound.append((pam_a,pam_b))
            bound = tuple(bound)
            pso_s = pso_scan(f, m_ffx, x_pam_start, var_name, bound,pam,f_x_flag ,pso_start_x_fx_dir, num_pso, POP_SIZE)
    #####################################################################################
    pso_start_x_fx_dir = {}
    f_x_flag = {}

    #
    xkpambest = os.getcwd() + os.sep + outdir + os.sep + 'xkpambest'
    if len(pam_) == 0:
        #
        xstart = ff('ffield', ['start'])
        xwrite = open(xkpambest, 'w')
        xwrite.write('')
        xwrite.close()
        pam = {}
        if scale_flag == 0:
            w_params(dir_path + os.sep + "ffield")
            paramename = 'params'
        else:
            paramename = 'params.scale'
        openparam = open(dir_path + os.sep + paramename, 'r')  #
        pamlines = openparam.readlines()
        openparam.close()
        numpam = 0
        for i in pamlines:
            perpam = i.split()
            if len(perpam) > 0 and i.find("#") != 0:
                numpam += 1
                pam[numpam] = {}
                pam[numpam]['id'] = {}
                pam[numpam]['value'] = {}
                pam[numpam]['id'][1] = int(perpam[0])
                pam[numpam]['id'][2] = int(perpam[1])
                pam[numpam]['id'][3] = int(perpam[2])
                pam[numpam]['value'][1] = float(perpam[-2])
                pam[numpam]['value'][2] = float(perpam[-1])
                pam[numpam]['step'] = float(perpam[-3])
    else:
        pam = pam_
    #
    for i in range(num_shunxu):
        pam_pso_scan(pam,f_x_flag ,pso_x_num, shunxu = 0)
        w_params(dir_path + os.sep + outdir + os.sep + "ffield.reax.start")
    for i in range(num_luanxu):
        pam_pso_scan(pam, f_x_flag, pso_x_num , shunxu = 1)
        w_params(dir_path + os.sep + outdir + os.sep + "ffield.reax.start")

if __name__ == "__main__":
    from ybs_ffield import ffield as ff
    f=ff
    psocan(ff, f, scale=0, e1=0.0001, E1=50.0, mul=0, outfile='outdata')
