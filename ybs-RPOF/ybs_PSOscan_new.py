import os
import math, re
import numpy as np
import matplotlib.pyplot as plt
from threading import Thread
import pandas as pd
from time import sleep
from ybs_ffieldpam import ffpam
from ybs_write_params import w_params
from ybs_mul_ffield import mul_ffield as m_ffx
import scipy.stats as stats
from ybs_id_v_flag_turn import id_v_flag_turn as pigturn
from ybs_pso_new import PSO

def PSOscan(ff,f,scale=0 ,n_dim = 5,size_pop = 10, max_iter = 10, num_shunxu = 1 , num_luanxu = 20, w = 0.8, c1 = 0.5, c2 = 0.5 ,pam_={},outfile='outdata'):
    #
    n_dim = n_dim
    size_pop = size_pop
    max_iter = max_iter
    scale_flag = scale
    outdir = outfile
    dir_path = os.getcwd()
    canshu_shezhi = [n_dim,size_pop,max_iter,w,c1,c2]

    ######
    def func(x_flag_value,ffield_name = ""):
        x_name = m_ffx(x_flag_value, ffield_name)
        f_v = f(x_name)[0]
        return f_v
    ############
    def pam_PSO_scan(pam,pam_step_list):
        n_dim = canshu_shezhi[0]
        size_pop = canshu_shezhi[1]
        max_iter = canshu_shezhi[2]
        can_w = canshu_shezhi[3]
        can_c1 = canshu_shezhi[4]
        can_c2 = canshu_shezhi[5]
        for pam_step in pam_step_list: #
            xx0_list = []
            a_list = []
            b_list = []
            x_id_list = []
            x_flag_value = []
            n_dim = len(pam_step)
            for num in pam_step:
                pam_a = pam[num]['value'][1]
                pam_b = pam[num]['value'][2]
                id1 = int(pam[num]['id'][1])
                id2 = int(pam[num]['id'][2])
                id3 = int(pam[num]['id'][3])
                pam_v = ffpam(xstart, [id1, id2, id3])
                #
                xx0_list.append(pam_v)
                a_list.append(pam_a)
                b_list.append(pam_b)
                x_id_list.append([id1,id2,id3])
                x_flag_value.append([id1,id2,id3,pam_v])
            popution = {}
            fx_name = ["PSO_new",]
            def mul_x_fun(xx): #
                ffield_name = fx_name[0]
                n_dim = len(pam_step)

                for i in range(n_dim):
                    xx_temp = xx[i]
                    x_flag_value[i][-1] =xx_temp
                xx_chrom = tuple(xx)
                if xx_chrom in popution:
                    f_v = popution[xx_chrom]
                else:
                    f_v = func(x_flag_value=x_flag_value,ffield_name=ffield_name)
                    popution[xx_chrom] = f_v
                return f_v
            #############################
            pso = PSO(func=mul_x_fun,dim=n_dim, pop=size_pop, max_iter= max_iter, X_x=xx0_list, lb=a_list, ub=b_list, w=can_w, c1=can_c1, c2=can_c2)
            pso.run()
            best_x = pso.gbest_x
            best_y = pso.gbest_y

            all_Y_history  =  pso.all_history_Y
            with open("Y_history",'w') as wy:
                for  y_history  in  all_Y_history:
                    y_history_list  = list(y_history[:,0].astype(str))

                    line  =  " ".join(y_history_list)
                    wy.write(line  +  "\n")

            '''
            plt.plot(pso.gbest_y_hist)
            plt.show()
            '''
            x_re_value = []
            for i in range(n_dim):
                pam_a = a_list[i]
                pam_b = b_list[i]
                id1 = x_id_list[i][0]
                id2 = x_id_list[i][1]
                id3 = x_id_list[i][2]
                pam_v = float("%.4f" % float(best_x[i]))
                x_re_value.append([id1,id2,id3,pam_v])
                with open(os.getcwd() + os.sep + outdir + os.sep + 'xkpambest', 'a') as w_pambest:
                    w_pambest.write(str(id1).ljust(3) + str(id2).ljust(3) + str(id3).ljust(5) \
                                    + str(pam_v).rjust(10) + str(pam_a).rjust(10) + str(pam_b).rjust(10) + '\n')
            with open(os.getcwd() + os.sep + outdir + os.sep + 'xkpambest', 'a') as w_pambest:
                w_pambest.write("#".ljust(50) + str(best_y) +'\n')
            pso_xf_best = m_ffx(x_re_value, "ffield.reax.start")
    #####################################################################################
    xkpambest = os.getcwd() + os.sep + outdir + os.sep + 'xkpambest'
    if len(pam_) == 0: #
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
    numpam_list = list(range(1,len(pam) + 1))
    pam_sunxu_step_list = [numpam_list[i:i+ n_dim] for i in range(0,len(numpam_list),n_dim)]
    #
    numpam_list2 = np.random.permutation(numpam_list)
    pam_suiji_step_list = [numpam_list2[i:i + n_dim] for i in range(0, len(numpam_list2), n_dim)]
    #
    # pso
    for i in range(num_shunxu):
        pam_PSO_scan(pam,pam_sunxu_step_list)
        w_params(dir_path + os.sep + outdir + os.sep + "ffield.reax.start")
    for i in range(num_luanxu):
        pam_PSO_scan(pam,pam_suiji_step_list)
        w_params(dir_path + os.sep + outdir + os.sep + "ffield.reax.start")
if __name__ == "__main__":
    pass



