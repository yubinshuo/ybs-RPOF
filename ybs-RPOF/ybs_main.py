# -*- coding:utf-8 -*-
import os, sys, time
import string
import math
import numpy as np
import matplotlib.pyplot as plt
from threading import Thread
from time import sleep
import sys
from ybs_bgf import bgfdata as wd
from ybs_trainset import trainset as ts
from ybs_lammpsderta import fderta as derta
from ybs_objfunction import Fx
from ybs_ffield import ffield as ff
from ybs_singlescan import sigscan
from ybs_ffieldpam import ffpam

def single_main(fcc=2,fgg=2,fee=2,fll=2,fhh=2,fmm = 0,llg=0):
    fc = fcc
    fg = fgg
    fe = fee
    fl = fll
    fh = fhh
    fm = fmm
    lg = llg

    tset, charge, geometry, energy,reaction_energy, cellpam, heat, nomini, mini = ts('trainset')
    geonatom, numdata, datatype, restraint, cell, geontype= wd('geofile')
    #optgeos, egeolist = optg()

    def fderta(x,lg,m):
        return derta(x,tset,datatype,restraint,nomini,mini,geonatom,cell,geontype,lg,m)

    def F(x,w1,w2,w3,w4,w5,w6,lg,m):
        return Fx(x,fderta(x,lg,m),tset,w1,w2,w3,w4,w5,w6)

    time_now = time.time()
    run_time = []
    run_time.append(time_now)
    def f(x,m=0):
        #x,wcharge,wgeo,wenergy cell heat
        start_time = time.time()
        run_time.append(start_time)
        run_sum_sec = int(run_time[-1]-run_time[0])
        run_min = run_sum_sec // 60
        run_sec = run_sum_sec % 60
        run_hour = run_min // 60
        if run_hour > 0:
            run_min = run_min % 60
        print("当前运行时间：%d hour %d min %d s" % (run_hour,run_min,run_sec))
        if len(run_time) > 2:
            print("当前 fx 计算运行时间：%.2f s" % (run_time[-1] - run_time[-2]))
        return F(x,fc,fg,fe,fl,fh,fm,lg,m)

    def dscan(e,E):
        return sigscan(ff,f,e,E)
    #####################################################################################
    #
    param_err_num = 0
    with open("params") as r_param:
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
                if pam_v < pam_a or pam_v > pam_b:
                    param_err_num += 1
                    print("pamrams file error: param %d-%d-%d value:%f out of range[%f %f] \n" %(pam_id1 ,pam_id2, pam_id3,pam_v,pam_a,pam_b) * 10)
    if param_err_num > 0:
        sys.exit()
    if not os.path.exists(os.curdir+os.sep+"outdata"):
        os.mkdir(os.curdir+os.sep+"outdata")

    xstart  = ff('ffield',['start'])
    fxstart = f(xstart)

    starterr = open(os.getcwd() +os.sep+"outdata" + os.sep + 'start.err.initial', 'w')
    starterr.close()
    starterr = open(os.getcwd() +os.sep+"outdata" +  os.sep + 'start.err.initial', 'a')
    starterr.write(" num Weight	Structure	QM	ForceField	Error \n")
    for key,value in fxstart[-1].items():
        starterr.write(str(key).center(7) +':' + value + '\n')
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

    diedai =dscan(0.0001,5.0)

    # 写best error
    dir_path = os.getcwd()
    dirlist = os.listdir(dir_path + os.sep + "outdata")
    fbestlist = []
    for i in dirlist:
        if i.startswith('ffield.reax.best'):
            fbestlist.append(i)
    for i in fbestlist:
        bestname = i[12:]
        print((bestname + "\n") * 100)
        fbest = f(("ffield.reax.%s" % bestname))
        besterr = open(os.getcwd() + os.sep + "outdata" + os.sep + 'start.err.%s' % bestname, 'w')
        besterr.close()
        besterr = open(os.getcwd() + os.sep + "outdata" + os.sep + 'start.err.%s' % bestname, 'a')
        besterr.write(
            " num Weight	Structure	                        QM	  FF.best    Error   FF.initial   Error\n")
        for key, value in fbest[6].items():
            num = key
            besterr.write(str(key).center(7) + ':' + value + str(fxstart[6][num].split()[-2]).rjust(10) + str(
                fxstart[6][num].split()[-1]).rjust(10) + '\n')
            if key == charge['num']:
                besterr.write("total error(Charge:initial,best) = %.4f, %.4f" % (fxstart[1], fbest[1]) + '\n')
            if key == (charge['num'] + geometry['num']):
                besterr.write("total error(Geometry:initial,best) = %.4f, %.4f" % (fxstart[2], fbest[2]) + '\n')
            if key == (charge['num'] + geometry['num'] + energy['num']):
                besterr.write("total error(Energy:initial,best) = %.4f, %.4f" % (fxstart[3], fbest[3]) + '\n')
            if key == (charge['num'] + geometry['num'] + energy['num'] + cellpam['num']):
                besterr.write("total error(Cellpam:initial,best) = %.4f, %.4f" % (fxstart[4], fbest[4]) + '\n')
            if key == (charge['num'] + geometry['num'] + energy['num'] + cellpam['num'] + heat['num']):
                besterr.write("total error(Heat:initial,best) = %.4f, %.4f" % (fxstart[5], fbest[5]) + '\n')
        besterr.write("TOTAL error:initial,best) = %.4f, %.4f" % (fxstart[0], fbest[0]) + '\n')
        besterr.close()
        # 作图
        readerror = open(os.getcwd() + os.sep + "outdata" + os.sep + 'start.err.%s' % bestname, 'r')
        rerror = readerror.readlines()

        for i in rerror:
            line = i
            if line.startswith("total error(Energy"):
                energyide = rerror.index(line)
                energyids = energyide - energy['num']
        errorlist = {}
        for i in rerror[energyids:energyide]:
            line = i
            spline = line.split()
            if 11 == len(spline):
                name1 = spline[3].split('/')[0]
                spname = name1.split('_')
                if len(spname) > 1:
                    if len(spname) < 3:
                        name = spname[0]
                    else:
                        name = spname[0] + '_' + spname[1]
                    if spname[-1] == "top":
                        value = int(rerror[rerror.index(line) - 1].split()[3].split('/')[0].split('_')[-1]) + 1
                    else:
                        value = int(spname[-1])
                    errorlist[name] = {}
        for i in rerror[energyids:energyide]:
            line = i
            spline = line.split()
            if 11 == len(spline):
                name1 = spline[3].split('/')[0]
                spname = name1.split('_')
                if len(spname) > 1:
                    if len(spname) < 3:
                        name = spname[0]
                    else:
                        name = spname[0] + '_' + spname[1]
                    if spname[-1] == "top":
                        value = int(rerror[rerror.index(line) - 1].split()[3].split('/')[0].split('_')[-1]) + 1
                    else:
                        value = int(spname[-1])
                    errorlist[name][value] = {}
        for i in rerror[energyids:energyide]:
            line = i
            spline = line.split()
            if 11 == len(spline):
                name1 = spline[3].split('/')[0]
                spname = name1.split('_')
                if len(spname) > 1:
                    if len(spname) < 3:
                        name = spname[0]
                    else:
                        name = spname[0] + '_' + spname[1]
                    if spname[-1] == "top":
                        value = int(rerror[rerror.index(line) - 1].split()[3].split('/')[0].split('_')[-1]) + 1
                    else:
                        value = int(spname[-1])
                    errorlist[name][value]['QM'] = float(spline[6])
                    errorlist[name][value]['FFstart'] = float(spline[-2])
                    errorlist[name][value]['FFbest'] = float(spline[-4])
        for k, v in errorlist.items():
            name = k
            QM = []
            FFbest = []
            FFstart = []
            x = []
            yname = 'Relative Energy (Kcal/mol)'
            if '_1a' in name or '_2a' in name or '_3a' in name or '_4a' in name or '_5a' in name or '_6a' in name or '_a' in name or \
                    '_7a' in name or '_8a' in name or '_9a' in name or '_10a' in name or '_11a' in name or '_12a' in name or \
                    '_13a' in name or '_14a' in name or '_15a' in name or '_16a' in name or '_17a' in name or '_18a' in name:
                xname = name + ' Valence Angle (degree)'
            elif '_b' in name or '_1b' in name or '_2b' in name or '_3b' in name or '_4b' in name or '_5b' in name or '_6b' in name or \
                    '_7b' in name or '_8b' in name or '_9b' in name or '_10b' in name or '_11b' in name or '_12b' in name or \
                    '_13b' in name or '_14b' in name or '_15b' in name or '_16b' in name or '_17b' in name or '_18b' in name:
                xname = name + ' Bond Distance (Å)'
            elif '_d' in name or '_1d' in name or '_2d' in name or '_3d' in name or '_4d' in name or '_5d' in name or '_6d' in name or \
                    '_7d' in name or '_8d' in name or '_9d' in name or '_10d' in name or '_11d' in name or '_12d' in name or \
                    '_13d' in name or '_14d' in name or '_15d' in name or '_16d' in name or '_17d' in name or '_18d' in name:
                xname = name + ' Dihedral Angle (degree)'
            elif '-v' in name:
                xname = name + ' Volume (100%)'
            elif '-1p' in name or '-2p' in name:
                xname = name + ' Times (Ps)'
            else:
                xname = name + ' Transtate'
            #
            if not os.path.exists(os.curdir + os.sep + "outdata" + os.sep + '%s' % bestname):
                os.mkdir(os.curdir + os.sep + "outdata" + os.sep + '%s' % bestname)
            bestpath = os.curdir + os.sep + "outdata" + os.sep + '%s' % bestname
            bngname = open(bestpath + os.sep + name, 'w')
            bngname.write('x'.rjust(10) + 'QM'.rjust(20) + 'FFstart'.rjust(20) + 'FFbest'.rjust(20) + '\n')
            bngname.close()
            for k1, v1 in v.items():
                if '_1a' in name or '_2a' in name or '_3a' in name or '_4a' in name or '_5a' in name or '_6a' in name or '_a' in name or \
                        '_7a' in name or '_8a' in name or '_9a' in name or '_10a' in name or '_11a' in name or '_12a' in name or \
                        '_13a' in name or '_14a' in name or '_15a' in name or '_16a' in name or '_17a' in name or '_18a' in name:
                    x.append(k1 * 1)
                elif '_b' in name or '_1b' in name or '_2b' in name or '_3b' in name or '_4b' in name or '_5b' in name or '_6b' in name or \
                        '_7b' in name or '_8b' in name or '_9b' in name or '_10b' in name or '_11b' in name or '_12b' in name or \
                        '_13b' in name or '_14b' in name or '_15b' in name or '_16b' in name or '_17b' in name or '_18b' in name:
                    x.append(k1 * 0.01)
                elif '_d' in name or '_1d' in name or '_2d' in name or '_3d' in name or '_4d' in name or '_5d' in name or '_6d' in name or \
                        '_7d' in name or '_8d' in name or '_9d' in name or '_10d' in name or '_11d' in name or '_12d' in name or \
                        '_13d' in name or '_14d' in name or '_15d' in name or '_16d' in name or '_17d' in name or '_18d' in name:
                    x.append(k1 * 1)
                elif '-v' in name:
                    x.append(k1 * 1)
                elif '-1p' in name:
                    x.append(k1 * 0.1)
                elif '-2p' in name:
                    x.append(k1 * 0.01)
                else:
                    x.append(k1)
                FFbest.append(v1['FFbest'])
                FFstart.append(v1['FFstart'])
                QM.append(v1['QM'])
                bngname = open(bestpath + os.sep + name, 'a')
                bngname.write(
                    str(k1).rjust(10) + str(v1['QM']).rjust(20) + str(v1['FFstart']).rjust(20) + str(
                        v1['FFbest']).rjust(
                        20) + '\n')
            bngname.close()
            sj = []
            sj = QM + FFbest + FFstart
            minqm = min(sj)
            maxqm = max(sj)
            mid = (minqm + maxqm) / 2
            if QM[-1] < mid or FFstart[-1] < mid or FFbest[-1] < mid:
                lengendid = "upper right"
            else:
                lengendid = "lower right"
            plt.plot(x, QM, 'o', color='k', linewidth=2.5, linestyle="-", label="QM")
            plt.plot(x, FFstart, 's', color='b', linewidth=2.5, linestyle="-", label="FFstart")
            plt.plot(x, FFbest, 'v', color='r', linewidth=2.5, linestyle="-", label="FFbest")
            # plt.title(name)
            plt.xlabel(xname)
            plt.ylabel(yname)
            plt.grid(True, linestyle='-.')
            plt.legend(loc=lengendid)
            plt.savefig(bestpath + os.sep + '%s.tif' % name, dpi=250, bbox_inches='tight')
            plt.close()
    print("恭喜程序运行完毕\n" * 10)
