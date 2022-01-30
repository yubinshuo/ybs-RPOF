# -*- coding:utf-8 -*-
import os, sys, time
import string
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib as mpl
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
def main_plot(fcc=2, fgg=2, fee=2, fll=2, fhh=2, fmm=2, llg=0,m=0):
    fc = fcc
    fg = fgg
    fe = fee
    fl = fll
    fh = fhh
    fm = fmm
    lg = llg
    m = m

    tset, charge, geometry, energy,reaction_energy, cellpam, heat, nomini, mini = ts('trainset')
    geonatom, numdata, datatype, restraint, cell, geontype = wd('geofile')

    def fderta(x, lg, m):
        return derta(x, tset, datatype, restraint, nomini, mini, geonatom, cell, geontype, lg, m)

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
        return F(x, fc, fg, fe, fl, fh, fm, lg, m)

    if not os.path.exists(os.curdir+os.sep+"outdata"):
        os.mkdir(os.curdir+os.sep+"outdata")
    lg_test_v = ffpam('ffield',[2,1,34])
    if lg_test_v == 999999 :
        lg = 0
    else:
        lg = 1

    xstart  = ff('ffield',['aaaa'])  #
    fxstart = f(xstart,m)

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
    dir_path = os.getcwd()
    dirlist     = os.listdir(dir_path+os.sep+"outdata" )
    fbestlist = []
    for i in dirlist:
        if i.startswith('ffield.reax.best'):
            fbestlist.append(i)
    for i in fbestlist:
        bestname = i[12:]
        print((bestname+"\n"))
        #
        lg_test_v = ffpam(("ffield.reax.%s"%bestname), [2, 1, 34])
        if lg_test_v == 999999:
            lg = 0
        else:
            lg = 1
        fbest = f(("ffield.reax.%s"%bestname),m)
        besterr = open(os.getcwd() +os.sep+"outdata" +  os.sep + 'start.err.%s'% bestname, 'w')
        besterr.close()
        besterr = open(os.getcwd() + os.sep+"outdata" + os.sep + 'start.err.%s'% bestname, 'a')
        besterr.write(" num Weight	Structure	                        QM	  FF.best    Error   FF.initial   Error\n")
        for key, value in fbest[-1].items():
            num = key
            besterr.write(str(key).center(7) + ':' + value + str(fxstart[-1][num].split()[-2]).rjust(10) + str(
                fxstart[-1][num].split()[-1]).rjust(10) + '\n')
            if key == charge['num']:
                besterr.write("total error(Charge:initial,best) = %.4f, %.4f" % (fxstart[1], fbest[1]) + '\n')
            if key == (charge['num'] + geometry['num']):
                besterr.write("total error(Geometry:initial,best) = %.4f, %.4f" % (fxstart[2], fbest[2]) + '\n')
            if key == (charge['num'] + geometry['num'] + energy['num']):
                besterr.write("total error(Energy:initial,best) = %.4f, %.4f" % (fxstart[3], fbest[3]) + '\n')
            # reaction_energy
            if key == (charge['num'] + geometry['num'] + energy['num'] + reaction_energy['num']):
                besterr.write("total error(Reaction_energy:initial,best) = %.4f, %.4f" % (fxstart[4], fbest[4]) + '\n')
            #
            if key == (charge['num'] + geometry['num'] + energy['num'] + reaction_energy['num'] + cellpam['num']):
                besterr.write("total error(Cellpam:initial,best) = %.4f, %.4f" % (fxstart[5], fbest[5]) + '\n')
            if key == (charge['num'] + geometry['num'] + energy['num'] + reaction_energy['num'] + cellpam['num'] + heat[
                'num']):
                besterr.write("total error(Heat:initial,best) = %.4f, %.4f" % (fxstart[6], fbest[6]) + '\n')
        besterr.write("TOTAL error:initial,best) = %.4f, %.4f" % (fxstart[0], fbest[0]) + '\n')
        besterr.close()
        zuobiao_biaoqian = {}
        try:
            with open(os.getcwd() + os.sep + "坐标标签") as r:
                zuobiaolines = r.readlines()
                for line in zuobiaolines:
                    if ":" in line:
                        spline = line.split(":")
                        name = spline[0]
                        biaoqian = spline[1]
                        zuobiao_biaoqian[name] = biaoqian
        except:
            pass
        #
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
                        value = int(rerror[rerror.index(line)-1].split()[3].split('/')[0].split('_')[-1]) + 1
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
                        value = int(rerror[rerror.index(line)-1].split()[3].split('/')[0].split('_')[-1]) + 1
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
                        value = int(rerror[rerror.index(line)-1].split()[3].split('/')[0].split('_')[-1]) + 1
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
            yname = 'Relative Energy (kcal/mol)'
            if '_1a' in name or '_2a' in name or '_3a' in name or '_4a' in name or '_5a' in name or '_6a' in name or '_a' in name or \
                    '_7a' in name or '_8a' in name or '_9a' in name or '_10a' in name or '_11a' in name or '_12a' in name or \
                    '_13a' in name or '_14a' in name or '_15a' in name or '_16a' in name or '_17a' in name or '_18a' in name:
                xname = name + '(valence scan)  angle (deg)'
            elif '_b' in name or '_1b' in name or '_2b' in name or '_3b' in name or '_4b' in name or '_5b' in name or '_6b' in name or \
                    '_7b' in name or '_8b' in name or '_9b' in name or '_10b' in name or '_11b' in name or '_12b' in name or \
                    '_13b' in name or '_14b' in name or '_15b' in name or '_16b' in name or '_17b' in name or '_18b' in name:
                xname =  name + '(bond scan)  Distance (Å)'
            elif '_d' in name or '_1d' in name or '_2d' in name or '_3d' in name or '_4d' in name or '_5d' in name or '_6d' in name or \
                    '_7d' in name or '_8d' in name or '_9d' in name or '_10d' in name or '_11d' in name or '_12d' in name or \
                    '_13d' in name or '_14d' in name or '_15d' in name or '_16d' in name or '_17d' in name or '_18d' in name:
                xname =  name + '(dihedral scan)  angle (deg)'
            elif '-v' in name :
                xname =  name + ' volume (%)'
            elif '-1p' in name or '-2p' in name:
                xname = name + ' Times (Ps)'
            elif '-m' in name or '-1m' in name or '-2m' in name or '-3m' in name or '-4m' in name or '-5m' in name:
                xname = name + ' Geometry scan in Cell'
            else:
                xname = name + '(bimolecular distance scan)  distance (Å)'
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
                elif '-v' in name :
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
                    str(k1).rjust(10) + str(v1['QM']).rjust(20) + str(v1['FFstart']).rjust(20) + str(v1['FFbest']).rjust(
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
            #
            #
            bwith = 1.0  #
            TK = plt.gca()  #
            TK.spines['bottom'].set_linewidth(bwith)  #
            TK.spines['left'].set_linewidth(bwith)  #
            TK.spines['top'].set_linewidth(bwith)  #
            TK.spines['right'].set_linewidth(bwith)  #
            ##
            font1 = {'family': 'Times New Roman',
                     'weight': 'normal',
                     'size': 12,
                     }
            font2 = {'family': 'Times New Roman',
                     'weight': 'normal',
                     'size': 12,
                     }

            #
            plt.plot(x, QM, 'o', color='k', linewidth=1.5,markersize=3, linestyle="-", label="DFT")
            plt.plot(x, FFstart, 's', color='b', linewidth=1.5,markersize=3, linestyle="-", label="ReaxFF" + "$_i$" + "$_n$" + "$_i$" + "$_t$" + "$_i$" + "$_a$" + "$_l$")
            plt.plot(x, FFbest, 'v', color='r', linewidth=1.5,markersize=3, linestyle="-", label="ReaxFF" + "$_N$" + "$_T$" + "$_O$")

            leg = plt.legend(labelspacing = 0.1,handletextpad=0.2,prop =font1)
            leg.get_frame().set_edgecolor('k')
            leg.get_frame().set_linewidth(0.5)

            if name in zuobiao_biaoqian:
                xname = zuobiao_biaoqian[name]
            #
            plt.xlabel(xname, font2)
            plt.ylabel(yname,font2)
            #
            plt.xticks(fontsize = 10)
            plt.yticks(fontsize=10)
            plt.grid(True, linestyle='-.')

            #
            plt.rcParams['figure.figsize'] = (5.0, 4.0)
            plt.rcParams['savefig.dpi'] = 400
            #ax.axis["right"].set_visible(False)
            plt.rc('font', family='Times New Roman')
            #
            plt.savefig(bestpath + os.sep + '%s.pdf' % name, bbox_inches='tight')
            plt.close()
