# -*- coding:utf-8 -*-
import os,sys
import time
from time import sleep
import timeout_decorator
from ybs_get_argv import get_argv
import string
import math
import numpy as np
import matplotlib.pyplot as plt
from lammps import lammps
import multiprocessing
import threading
from ybs_bgf import bgfdata as wd
from ybs_trainset import trainset as ts
from ybs_mpiderta import fderta as derta
from ybs_objfunction import Fx
from ybs_ffield import ffield as ff
from ybs_ffieldpam import ffpam
from ybs_lmp_run import ybs_lmpnomini,ybs_lmpmini
from ybs_singlescan_func import sigscan_func as sigfunc
from ybs_singlescan import sigscan
from ybs_SAscan import SAscan
from ybs_psoscan import psoscan as psoscan
from ybs_GAscan import GAscan
from ybs_PSOscan_new import PSOscan as PSOscan_new
from ybs_reassign_geolist import reassign_geo

outdir = "outdata"
argv = sys.argv
os.makedirs(os.path.join(os.getcwd(),outdir), exist_ok=True)

## get argv
get_argv_return = get_argv(argv)
mpi_nc, mnp, md, m, scale_flag = get_argv_return[0]
single_scan_flag, mul_pso, pso_flag, pso_new_flag, sa_flag, ga_flag, paowu_plus_pso_flag, mini_all_flag = get_argv_return[1]
fc, fg, fe, fl, fh, fm = get_argv_return[2]
fix_error = get_argv_return[3]
##
lg = 0
num_shunxu, num_luanxu = 1, 20
n_dim, size_pop, precision, max_iter = 10, 6, 1e-4, 100
# pso
pso_c1, pso_c2, pso_w = 1.5, 1.5, 0.2
# Sa
T_max, T_min, sa_L, max_stay_counter = 1, 1e-6, 2, 5
pam = {}
num_pso = 1
#
if mpi_nc == 1 and "-p" not in argv:
    from ybs_main import single_main
    main = single_main(fc,fg,fe,fl,fh,fm,lg)
elif "-p" in argv:
    from  ybs_mainplot import main_plot
    m_plot = main_plot(fc,fg,fe,fl,fh,fm,lg,mini_all_flag)
else:
    tset, charge, geometry, energy, reaction_energy, cellpam, heat, nomini, mini = ts('trainset')
    geonatom, numdata, datatype, restraint, cell, geontype= wd('geofile')
    nomini_single = nomini.copy()
    mini_single  = mini.copy()
    [chirld_mini,chirld_nomini,chirld_nomini_cell,mpi_n,mpi_m] = reassign_geo(mpi_nc, mnp, mini, nomini, cell, md)

    def lmpnomini(ffi,nom,m=0):
        return ybs_lmpnomini(ffi,nom,restraint,geonatom,datatype,cell,lg,m,outdir='outdata')
    def lmpmini(ffi, nom, m=0):
        return ybs_lmpmini(ffi,nom,restraint,geonatom,datatype,cell,lg,m,outdir='outdata')
    ### define mpirun
    def mpi_nomini(ff,m=0):
        ffield = ff
        mi = m
        result = []
        fff = {}
        threads = []
        if len(nomini) != 0:
            for i in range(mpi_n):
                threads.append(threading.Thread(target =lmpnomini, args=(ffield,chirld_nomini[i],mi,)))
            if __name__ == '__main__':
                pool = multiprocessing.Pool(processes=mpi_n)
                for i in range(mpi_n):
                    result.append(pool.apply_async(lmpnomini,(ffield,chirld_nomini[i],mi,)))
                pool.close()
                pool.join()
                #for t in threads:
                    #t.setDaemon(True)
                    #t.start()
                #t.join()
            for res in result:
                re_dict = dict(res.get())
                for kk,vv in re_dict.items():
                    fff[kk] = vv
        return fff
    ### define mpirun mini
    def mpi_mini(ff,m=0):
        ffield = ff
        mi = m
        result = []
        fff = {}
        threads = []
        if len(mini) != 0:
            for i in range(mpi_m):
                threads.append(threading.Thread(target =lmpmini, args=(ffield,chirld_mini[i],mi,)))
            if __name__ == '__main__':
                pool = multiprocessing.Pool(processes=mpi_m)
                for i in range(mpi_m):
                    result.append(pool.apply_async(lmpmini,(ffield,chirld_mini[i],mi,)))
                pool.close()
                pool.join()
                #for t in threads:
                    #t.setDaemon(True)
                    #t.start()
                #t.join()
            for res in result:
                re_dict = dict(res.get())
                for kk,vv in re_dict.items():
                    fff[kk] = vv
        return fff
    ## define mpi_sum
    def mpi_sum(ff, m=0):
        ffield = ff
        mi = m
        result = []
        fff1 = {}
        fff2 = {}
        threads = []
        for i in range(mpi_nc):
            mpi_num = i
            if mpi_num < mpi_n:
                threads.append(threading.Thread(target=lmpnomini, args=(ffield, chirld_nomini[i], mi,)))

            else:
                i2 = i - mpi_n
                threads.append(threading.Thread(target=lmpmini, args=(ffield, chirld_mini[i2], mi,)))

        if __name__ == '__main__':
            pool = multiprocessing.Pool(processes=mpi_nc)
            for i in range(mpi_nc):
                mpi_num = i
                if mpi_num < mpi_n:
                    result.append(pool.apply_async(lmpnomini, (ffield, chirld_nomini[i], mi,)))
                else:
                    i2 = i - mpi_n
                    result.append(pool.apply_async(lmpmini, (ffield, chirld_mini[i2], mi,)))

            pool.close()
            pool.join()
            # for t in threads:
            # t.setDaemon(True)
            # t.start()
            # t.join()
        for res in result:
            res_index = result.index(res)
            if res_index < mpi_n:
                re_dict1 = dict(res.get())

                for kk1, vv1 in re_dict1.items():
                    fff1[kk1] = vv1
            else:
                re_dict2 = dict(res.get())
                for kk2, vv2 in re_dict2.items():
                    fff2[kk2] = vv2
        return fff1, fff2
    def fderta(x,m):
        return derta(x,tset,geontype,mpi_sum(x,m))
    def F(x,w1,w2,w3,w4,w5,w6,m):
        return Fx(x,fderta(x,m),tset,w1,w2,w3,w4,w5,w6,fix_error)
    time_now = time.time()
    run_time = []
    run_time.append(time_now)
    num_fun = 0
    num_single_fun = 0
    def fxxx(x,m=0):
        #x,wcharge,wgeo,wenergy
        global num_fun
        global num_single_fun
        num_fun += 1
        start_time = time.time()
        run_time.append(start_time)
        run_sum_sec = int(run_time[-1]-run_time[0])
        run_min = run_sum_sec // 60
        run_sec = run_sum_sec % 60
        run_hour = run_min // 60
        if run_hour > 0:
            run_min = run_min % 60
        if len(run_time) > 2:
            pass
        return F(x,fc,fg,fe,fl,fh,fm,m)
    #--------------------------------------- single ---------------------------------------#
    def fderta_single(x,m):
        print("yes_single")
        return derta(x,tset,geontype,[lmpnomini(x, nomini_single,m),lmpmini(x, mini_single, m)])
    def F_single(x,w1,w2,w3,w4,w5,w6,m):
        return Fx(x,fderta_single(x,m),tset,w1,w2,w3,w4,w5,w6,fix_error)
    def f_single(x,m=0):
        #x,wcharge,wgeo,wenergy cell heat
        global num_fun
        global num_single_fun
        num_single_fun += 1
        start_time = time.time()
        run_time.append(start_time)
        run_sum_sec = int(run_time[-1]-run_time[0])
        run_min = run_sum_sec // 60
        run_sec = run_sum_sec % 60
        run_hour = run_min // 60
        if run_hour > 0:
            run_min = run_min % 60
        if len(run_time) > 2:
            pass
        return F_single(x,fc,fg,fe,fl,fh,fm,m)
    ######
    if "-pmpi" in argv:
        n_time_run = 999999999
    else:
        n_time_run = 200
    #
    @timeout_decorator.timeout(n_time_run)
    def f(x):
        try:
            #return f_single(x, m)
            print("m=",m)
            return fxxx(x, m)
        except:
            return f_single(x,m)
    def dscan(e,E):
        return sigscan(ff,f,scale_flag,e,E,mul_pso)
    #tset eofile
    geonatom, numdata, datatype, restraint, cell, geontype= wd('geofile')
    all_geo_in_tset = nomini + mini
    for  geo_name_test in  all_geo_in_tset:
        try:geontype[geo_name_test]
        except: print(geo_name_test," "),sys.exit()
    #
    param_err_num = 0
    with open("params.scale") as r_param:
        numpam = 0
        pam = {}
        pamlines = r_param.readlines()
        for i in pamlines:
            spam = i.split()
            if len(spam) > 0 and i.find("#") != 0:
                pam_id1 = int(spam[0])
                pam_id2 = int(spam[1])
                pam_id3 = int(spam[2])
                pam_a = float(spam[-2])
                pam_b = float(spam[-1])
                pam_v = ffpam("ffield", [pam_id1, pam_id2, pam_id3])
                numpam += 1
                pam[numpam] = {}
                pam[numpam]['id'] = {}
                pam[numpam]['value'] = {}
                pam[numpam]['id'][1] = int(spam[0])
                pam[numpam]['id'][2] = int(spam[1])
                pam[numpam]['id'][3] = int(spam[2])
                pam[numpam]['value'][1] = pam_a
                pam[numpam]['value'][2] = pam_b
                pam[numpam]['step'] = float(spam[-3])
                if pam_v < pam_a or pam_v > pam_b:
                    param_err_num += 1
                    print("pamrams file error: param %d-%d-%d value:%f out of range[%f %f] \n" %(pam_id1 ,pam_id2, pam_id3,pam_v,pam_a,pam_b) * 10)
    if param_err_num > 0:
        sys.exit()
    #
    lg_test_v = ffpam('ffield',[2,1,34])
    if lg_test_v == 999999 :
        lg = 0
    else:
        lg = 1
    xstart  = ff('ffield',['start'])
    fxstart = f(xstart)
    starterr = open(os.getcwd() +os.sep+"outdata" + os.sep + 'start.err.initial', 'w')
    starterr.close()
    starterr = open(os.getcwd() +os.sep+"outdata" +  os.sep + 'start.err.initial', 'a')
    starterr.write(" num Weight	Structure	QM	ForceField	Error \n")
    for key,value in fxstart[-1].items():
        starterr.write(str(key).center(7) +' : ' + value + '\n')
        if key == charge['num']:
            starterr.write("total error(Charge) = %.4f"% fxstart[1] + '\n' )
        if key == (charge['num']+geometry['num']):
            starterr.write("total error(Geometry) = %.4f" % fxstart[2] + '\n')
        if key == (charge['num']+geometry['num']+energy['num']):
            starterr.write("total error(Energy) = %.4f" % fxstart[3] + '\n')
        #reaction_energy
        if key == (charge['num'] + geometry['num'] + energy['num'] + reaction_energy['num']):
            starterr.write("total error(Reaction_energy = %.4f" % fxstart[4] + '\n')
        #
        if key == (charge['num']+geometry['num']+energy['num']+ reaction_energy['num']+cellpam['num']):
            starterr.write("total error(Cell) = %.4f" % fxstart[5] + '\n')
        if key == (charge['num']+geometry['num']+energy['num']+ reaction_energy['num']+cellpam['num']+heat['num']):
            starterr.write("total error(Heat) = %.4f" % fxstart[6] + '\n')
    starterr.write("TOTAL error = %.4f" % fxstart[0] + '\n')
    starterr.close()
    if "-pmpi" in argv:
        from ybs_plot_mpi import mpi_plot
        fbestlist = []
        dir_path = os.getcwd()
        dirlist = os.listdir(dir_path + os.sep + "outdata")
        for i in dirlist:
            if i.startswith('ffield.reax.best'):
                fbestlist.append(i)
        for i in fbestlist:
            bestname = i[12:]
            print((bestname + "\n"))
            lg_test_v = ffpam(("ffield.reax.%s" % bestname), [2, 1, 34])
            if lg_test_v == 999999:
                lg = 0
            else:
                lg = 1
            fbest = f(("ffield.reax.%s" % bestname))
            mpi_plot(fxstart, fbest, bestname)
    else:
        #
        if single_scan_flag == 1:
            diedai = dscan(0.0001,5.0)
        elif pso_flag == 1:
            with open(os.path.join(outdir,"pso_result"), 'w') as ww:
                ww.write("")
            psoscan(ff, f, scale_flag, num_shunxu, num_luanxu, pso_c1, pso_c2, pso_w, num_pso, size_pop, n_dim,pam)
        elif paowu_plus_pso_flag == 1:
            for i in range(100):
                sigfunc(ff, f, scale_flag, pam, K=1, mul=0)
                num_shunxu = 1
                num_luanxu = 1
                psoscan(ff, f, scale_flag, num_shunxu, num_luanxu, pso_c1, pso_c2, pso_w, num_pso, size_pop, n_dim, pam)
        elif ga_flag == 1:
            GAscan(ff,f,scale_flag,n_dim ,size_pop ,precision, max_iter,num_shunxu , num_luanxu)
        elif pso_new_flag == 1:
            PSOscan_new(ff, f, scale_flag, n_dim, size_pop, max_iter , num_shunxu, num_luanxu , pso_w, pso_c1, pso_c2)
        elif sa_flag == 1:
            SAscan(ff, f, scale_flag, n_dim , T_max, T_min, sa_L, max_stay_counter, num_shunxu, num_luanxu)