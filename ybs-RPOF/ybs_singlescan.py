import os
import math
import numpy as np
import matplotlib.pyplot as plt
from threading import Thread
from time import sleep
from ybs_ffieldpam import ffpam
from ybs_write_params import w_params
from ybs_mul_ffield import mul_ffield as m_ffx
from ybs_pso import PSO
def sigscan(ff,f,scale=0,e1=0.0001,E1=50.0,mul=0,outfile='outdata'):
    global K,xstart,e,E
    E=E1
    e=e1
    scale_flag = scale
    outdir = outfile
    dir_path = os.getcwd()
    mul_pso = mul
    pam = {}
    pam_start = {}
    pam_final = {}
    XX = {}
    YY = {}
    ### params.scale
    pam_scale = {}
    with open(dir_path + os.sep + "params.scale") as r:
        lines = r.readlines()
        for i in lines:
            perpam = i.split()
            if len(perpam) > 0 and i.find("#") != 0:
                pam_flag = perpam[0]+'-'+perpam[1]+'-'+perpam[2]
                pam_scale[pam_flag]={}
                pam_scale[pam_flag]['a'] = float(perpam[-2])
                pam_scale[pam_flag]['b'] = float(perpam[-1])

    ### params
    if scale_flag == 0:
        write_params = w_params(dir_path + os.sep + "ffield")
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
            XX[numpam] = {}
            YY[numpam] = {}
    Npam = numpam
    xstart = ff('ffield',['start'])

    EF = {}
    K = 0
    nn = 200
    fxstart = f(xstart)
    EF[0] = fxstart[0]

    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'w')
    xxxx.write('')
    xxxx.close()
    xwrite = open(os.getcwd() + os.sep + outdir + os.sep + 'xkpambest', 'w')
    xwrite.write('')
    xwrite.close()
    error = open(os.getcwd() + os.sep + outdir + os.sep + 'ferror', 'w')
    error.write(str(K).ljust(5) + str(EF[0]) + '\n')
    error.close()
    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + ' ', 'a')
    xxxx.write('')
    xxxx.close()
    dertaK = 2
    #
    nK = 0
    NK = {}
    NK[0] = 0
    #
    FFFF = 0
    def pso_scan(f, m_ffx, x_pam_start, var_name, bound,pso_start_x_fx_dir,num_pso=1,POP_SIZE=10):
        print("yes")
        f = f
        m_ffx = m_ffx
        x_pam_start = x_pam_start
        var_name = var_name
        bound = bound
        POP_SIZE = POP_SIZE
        w = 1
        c1 = 0.5
        c2 = 0.5
        v_max = bound[0][1]-bound[0][0]
        x_fx_list = {}
        f_x_flag = {}
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
        ############# def
        def func(*x_list):
            x_flag_value = {}
            x_flag_value_str_list = ""
            for i in range(len(var_name)):
                pam_flag = var_name[i]
                pam_v = float("%.4f" % x_list[i])
                x_flag_value[pam_flag] = pam_v
                print(pam_flag, pam_v)
                x_flag_value_str = pam_flag + "_" + str(pam_v) + "_"
                x_flag_value_str_list += x_flag_value_str
            pso_x_name = m_ffx(x_flag_value, "pso")
            if x_flag_value_str_list in f_x_flag:
                # x_flag_value_num += 1
                # print(x_flag_value_str_list + "___" + str(f_x_flag[x_flag_value_str_list]) +  "___" + str(x_flag_value_num))
                f_v = -f_x_flag[x_flag_value_str_list]
            else:
                f_v = -f(pso_x_name)[0]
                f_x_flag[x_flag_value_str_list] = -f_v
                # 
                x_fx_list[pam_v] = -f_v
                with open("outdata/pso_cal", 'a') as www:
                    www.write((x_flag_value_str_list + "  ").ljust(40) + str(-f_v) + "\n")
            return f_v
        ####################### end def
        pso = PSO(func, bound, POP_SIZE, x_pam_start, w, c1, c2, v_max)
        for _ in range(num_pso):
            pso.pso()
        ##
        with open("outdata" + os.sep + "pso_result") as r:
            lines = r.readlines()
            min_f_v = 99999999999999
            min_x_v = []
            x_re_value = {}
            for line in lines:
                spline = line.split()
                if len(spline) != 0:
                    pso_re_pam_val_list = spline[0].replace("[", '').replace("]", "").split(",")
                    pso_re_fx_v = float(spline[1])
                    if pso_re_fx_v < min_f_v:
                        min_f_v = pso_re_fx_v
                        min_x_v = pso_re_pam_val_list
        for i in range(len(var_name)):
            pam_flag = var_name[i]
            pam_v = float("%.4f" % float(min_x_v[i]))
            x_re_value[pam_flag] = pam_v
        # pso_x_name = m_ffx(x_re_value, "best")
        return x_re_value,x_fx_list ,min_f_v
    #####################################################################################
    while dertaK > 1:
        ##
        K += 1
        numpam = 0
        pam_start = pam
        tt = {}
        ss = {}
        if K in range(100):
            pam_ = pam_start
        else:
            pam_ = pam_final
        xwrite = open(os.getcwd() + os.sep + outdir + os.sep + 'xkpambest', 'a')
        xwrite.write('\n' + 'K=' + str(K).rjust(5) + '\n')
        xwrite.close()

        re_pam = {}
        re_num = 0
        for key,val in pam_.items():
            numpam = key
            k = 0
            a = pam[numpam]['value'][1]
            b = pam[numpam]['value'][2]
            X2, Y = list(), list()
            XX[numpam][K] = []
            YY[numpam][K] = []
            id1 = int(pam[numpam]['id'][1])
            id2 = int(pam[numpam]['id'][2])
            id3 = int(pam[numpam]['id'][3])
            pam_flag = str(pam[numpam]['id'][1])+'-'+str(pam[numpam]['id'][2])+'-'+str(pam[numpam]['id'][3])
            x2s = {}
            pso_start_x_fx_dir = {}
            #
            def xf(x,name):
                global xstart
                return ff(xstart,[id1,id2,id3,x,str(name)])
            #
            def f_pwx(x,pso_start_x_fx_dir):
                pwx_x_name = xf(x, "pwx")
                if x in pso_start_x_fx_dir:
                    f_v = pso_start_x_fx_dir[x]
                    print(" ")
                else:
                    f_v = f(pwx_x_name)[0]
                    pso_start_x_fx_dir[x] = f_v
                return f_v

            x_start = pam_scale[pam_flag]['a']
            x_end = pam_scale[pam_flag]['b']
            #
            x1 = float('%.4f'% a)
            x3 = float('%.4f'% b)
            #x2 = float('%.4f'% ((a + b) / 2.0))
            ee = e + 1
            xx1 = xf(x1, 'x1')
            #xx2 = xf(x2, 'x2')
            xx3 = xf(x3, 'x3')
            x_left = []
            x_right = []
            #x1
            xxxx= open(os.getcwd() +os.sep+outdir + os.sep + 'xxxx', 'a')
            xxxx.write('\n'+ str(K).ljust(3)+ str(numpam).ljust(5) +str(id1).ljust(3) + str(id2).ljust(3) + str(id3).ljust(3) + ' x1=%s' % (x1))
            xxxx.close()
            fx1 = float('%.6f'% f_pwx(x1,pso_start_x_fx_dir))
            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
            xxxx.write(' fx1=%s'%(fx1))
            xxxx.close()
            x_a = x1
            f_a = fx1
            x2s[x1] = fx1

            #x3
            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
            xxxx.write( ' x3=%s' % (x3))
            xxxx.close()
            fx3 = float('%.6f'% f_pwx(x3,pso_start_x_fx_dir))
            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
            xxxx.write(' fx3=%s' % (fx3))
            xxxx.close()
            x_b = x3
            f_b = fx3
            x2s[x3] = fx3

            #s=0.0 x2
            s0 = ffpam(xstart, [id1, id2, id3])
            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
            xxxx.write('\n' + ' x0=%s' % (s0))
            xxxx.close()
            xs0 = xf(s0, 'x2')
            fs0 = float('%.6f'% f_pwx(s0,pso_start_x_fx_dir))
            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
            xxxx.write( ' fx0=%s' % (fs0))
            xxxx.close()
            x2s[s0] = fs0
            x_left.append((x1,fx1))
            x_right.append((x3,fx3))

            x_pam_start = [s0]
            var_name = [pam_flag]
            bound = ((x_a,x_b),)
            #  pso
            pso_start_x_fx_dir.update(x2s)
            if mul_pso == 1:
                pso_s = pso_scan(f, m_ffx, x_pam_start, var_name, bound,pso_start_x_fx_dir,num_pso=1,POP_SIZE=8)
                pso_start_x_fx_dir.update(pso_s[1])

                x2 = pso_s[0][pam_flag]
                fx2 = pso_s[-1]
                pso_x_fx_list = sorted(pso_s[1].items(), key=lambda item: item[1])
                xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                xxxx.write('\n The start x2=%s' % (x2))
                xxxx.close()
                xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                xxxx.write(' fx2=%s                               PSO' % (fx2))
                xxxx.close()
                if x2 == x_start or x2 == x_end:
                    x2 = pso_x_fx_list[1][0]
                    fx2 = pso_x_fx_list[1][1]
                    pso_x_fx_list = sorted(pso_s[1].items(), key=lambda item: item[1])
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write('\n The start x2=%s' % (x2))
                    xxxx.close()
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write(' fx2=%s                               PSO' % (fx2))
                    xxxx.close()

            else:
                if fs0 > fx1 or fs0 > fx3:
                    pso_s = pso_scan(f, m_ffx, x_pam_start, var_name, bound, pso_start_x_fx_dir, num_pso=1, POP_SIZE=8)
                    pso_start_x_fx_dir.update(pso_s[1])

                    x2 = pso_s[0][pam_flag]
                    fx2 = pso_s[-1]
                    pso_x_fx_list = sorted(pso_s[1].items(), key=lambda item: item[1])
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write('\n The start x2=%s' % (x2))
                    xxxx.close()
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write(' fx2=%s                               PSO' % (fx2))
                    xxxx.close()
                    if x2 == x_start or x2 == x_end:
                        x2 = pso_x_fx_list[1][0]
                        fx2 = pso_x_fx_list[1][1]
                        pso_x_fx_list = sorted(pso_s[1].items(), key=lambda item: item[1])
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write('\n The start x2=%s' % (x2))
                        xxxx.close()
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write(' fx2=%s                               PSO' % (fx2))
                        xxxx.close()
                else:
                    x2 = s0
                    fx2 = fs0
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write('\n The start x2=%s' % (x2))
                    xxxx.close()
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write(' fx2=%s  ' % (fx2))
                    xxxx.close()
            if x2 == x_start or x2 == x_end:
                x2 = s0
                fx2 = fs0

            fx = fx2
            x = x2
            x_list = []
            x2m = {}
            while ee > e:
                k += 1
                if fx2 > fx1 or fx2 > fx3:
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + '凸点', 'a')
                    xxxx.write('\n' + ('%s %s %s %s %s %s %s %s %s' % (K,numpam,id1,id2,id3,k,fx1,fx2,fx3)))
                    xxxx.close()
                    fmmin = min(fx1,fx3)

                    #m0 = 0.0001
                    m0 = x1 + 0.0005
                    xm0 = xf(m0, 'x2')
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write('\n' + (' x2=%s' % m0))
                    xxxx.close()
                    fm0 = float('%.6f'% f_pwx(m0,pso_start_x_fx_dir))
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write((' fx2=%s\n' % fm0))
                    xxxx.close()
                    x2m[m0] = fm0

                    #m1 = 0.0125
                    m1 = float('%.4f'% (x1 + 0.0125 * (x3 - x1)))
                    xm1 = xf(m1, 'x2')
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write('\n' + (' x2=%s' % m1))
                    xxxx.close()
                    fm1 = float('%.6f'% f_pwx(m1,pso_start_x_fx_dir))
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write((' fx2=%s\n' % fm1))
                    xxxx.close()
                    x2m[m1] = fm1

                    #m3 = 0.05
                    m3 = float('%.4f'% (x1 + 0.05 * (x3 - x1)))
                    xm3 = xf(m3, 'x2')
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write('\n' + (' x2=%s' % m3))
                    xxxx.close()
                    fm3 = float('%.6f'% f_pwx(m3,pso_start_x_fx_dir))
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write((' fx2=%s\n' % fm3))
                    xxxx.close()
                    x2m[m3] = fm3

                    #m5 = 0.1
                    m5 = float('%.4f'% (x1 + 0.1 * (x3 - x1)))
                    xm5 = xf(m5, 'x2')
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write('\n' + (' x2=%s' % m5))
                    xxxx.close()
                    fm5 = float('%.6f'% f_pwx(m5,pso_start_x_fx_dir))
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write((' fx2=%s\n' % fm5))
                    xxxx.close()
                    x2m[m5] = fm5

                    #m7 = 0.9
                    m7 = float('%.4f'% (x1 + 0.9 * (x3 - x1)))
                    xm7 = xf(m7, 'x2')
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write('\n' + (' x2=%s' % m7))
                    xxxx.close()
                    fm7 = float('%.6f'% f_pwx(m7,pso_start_x_fx_dir))
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write((' fx2=%s\n' % fm7))
                    xxxx.close()
                    x2m[m7] = fm7

                    #m9 = 0.995
                    m9 = float('%.4f'% (x1 + 0.995 * (x3 - x1)))
                    xm9 = xf(m9, 'x2')
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write('\n' + (' x2=%s' % m9))
                    xxxx.close()
                    fm9 = float('%.6f'% f_pwx(m9,pso_start_x_fx_dir))
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write((' fx2=%s\n' % fm9))
                    xxxx.close()
                    x2m[m9] = fm9
                    pso_start_x_fx_dir.update(x2m)
                    x2mmin = sorted(x2m.items(),key=lambda item:item[1])[0][0]
                    x2fmmin = sorted(x2m.items(),key=lambda item:item[1])[0][1]
                    if x2fmmin < fmmin:
                        x2 = x2mmin
                        fx2 = x2fmmin
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write('\n' + (' x2=%s' % x2))
                        xxxx.close()
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write((' fx2=%s\n' % fx2))
                        xxxx.close()
                    else:
                        x = ffpam(xstart, [id1, id2, id3])
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write('\n' + (' x=%s' % x))
                        xxxx.close()
                        xx = xf(x, 'x')
                        fx = float('%.6f' % f_pwx(x,pso_start_x_fx_dir))
                        fx2 = fx
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write((' fx=%s\n' % fx))
                        xxxx.close()
                        breakx = xf(x, ('break_%s_%s-%s-%s_%s' % (numpam, id1, id2, id3, K)))
                        break
                elif k==1 and (fx1 == fx2 or fx2 ==fx3):
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + ' ', 'a')
                    xxxx.write('\n' + ('%s %s %s %s %s %s %s %s %s' % (K,numpam,id1,id2,id3,k,fx1,fx2,fx3)))
                    xxxx.close()
                    if fx1 == fx2 == fx3:
                        x = ffpam(xstart,[id1,id2,id3])
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write('\n' + (' x=%s' % x))
                        xxxx.close()
                        fx = fx1
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write((' fx=%s\n' % fx))
                        xxxx.close()
                        print("  ")
                        breakx3 = xf(x1, ('break_%s_%s-%s-%s_%s' % (numpam, id1, id2, id3, K)))
                        break
                    elif fx1 > fx2:
                        x = x2
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write('\n' + (' x=%s' % x))
                        xxxx.close()

                        fx = fx2
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write((' fx=%s\n' % fx))
                        xxxx.close()
                        breakx = xf(x, ('break_%s_%s-%s-%s_%s' % (numpam, id1, id2, id3, K)))
                        break
                    elif fx3 > fx2:
                        x = x2
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write('\n' + (' x=%s' % x))
                        xxxx.close()

                        fx = fx2
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write((' fx=%s\n' % fx))
                        xxxx.close()
                        breakx = xf(x, ('break_%s_%s-%s-%s_%s' % (numpam, id1, id2, id3, K)))
                        break
                    elif fx1 < fx2:
                        x3 = x2
                        fx3 = fx2
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write('\n' + (' x3=%s' % x3))
                        xxxx.close()
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write((' fx3=%s\n' % fx3))
                        xxxx.close()
                        x2m = {}

                        # m0 = 0.0001
                        m0 = x1 + 0.0005
                        xm0 = xf(m0, 'x2')
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write('\n' + (' x2=%s' % m0))
                        xxxx.close()
                        fm0 = float('%.6f'% f_pwx(m0,pso_start_x_fx_dir))
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write((' fx2=%s\n' % fm0))
                        xxxx.close()
                        x2m[m0] = fm0

                        # m1=0.1
                        m1 = float('%.4f'% (x1 + 0.1 * (x3 - x1)))
                        xm1 = xf(m1, 'x2')
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write('\n' + (' x2=%s' % m1))
                        xxxx.close()
                        fm1 = float('%.6f'% f_pwx(m1,pso_start_x_fx_dir))
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write((' fx2=%s\n' % fm1))
                        xxxx.close()
                        x2m[m1] = fm1

                        # m2=0.2
                        m2 = float('%.4f'% (x1 + 0.2 * (x3 - x1)))
                        xm2 = xf(m2, 'x2')
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write('\n' + (' x2=%s' % m2))
                        xxxx.close()
                        fm2 = float('%.6f'% f_pwx(m2,pso_start_x_fx_dir))
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write((' fx2=%s\n' % fm2))
                        xxxx.close()
                        x2m[m2] = fm2

                        # m3=0.3
                        m3 = float('%.4f'% (x1 + 0.3 * (x3 - x1)))
                        xm3 = xf(m3, 'x2')
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write('\n' + (' x2=%s' % m3))
                        xxxx.close()
                        fm3 = float('%.6f'% f_pwx(m3,pso_start_x_fx_dir))
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write((' fx2=%s\n' % fm3))
                        xxxx.close()
                        x2m[m3] = fm3

                        # m4 = 0.4
                        m4 = float('%.4f'% (x1 + 0.4 * (x3 - x1)))
                        xm4 = xf(m4, 'x2')
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write('\n' + (' x2=%s' % m4))
                        xxxx.close()
                        fm4 = float('%.6f'% f_pwx(m4,pso_start_x_fx_dir))
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write((' fx2=%s\n' % fm4))
                        xxxx.close()
                        x2m[m4] = fm4

                        # m5 = 0.5
                        m5 = float('%.4f'% (x1 + 0.5 * (x3 - x1)))
                        xm5 = xf(m5, 'x2')
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write('\n' + (' x2=%s' % m5))
                        xxxx.close()
                        fm5 = float('%.6f'% f_pwx(m5,pso_start_x_fx_dir))
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write((' fx2=%s\n' % fm5))
                        xxxx.close()
                        x2m[m5] = fm5
                        pso_start_x_fx_dir.update(x2m)
                        x2mmin = sorted(x2m.items(), key=lambda item: item[1])[0][0]
                        x2fmmin = sorted(x2m.items(), key=lambda item: item[1])[0][1]
                        if x2fmmin < fx1:
                            x2 = x2mmin
                            fx2 = x2fmmin
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write('\n' + (' x2=%s' % x2))
                            xxxx.close()
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write((' fx2=%s\n' % fx2))
                            xxxx.close()
                        else:
                            x = ffpam(xstart, [id1, id2, id3])
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write('\n' + (' x=%s' % x))
                            xxxx.close()
                            xx = xf(x, 'x')
                            fx = float('%.6f' % f_pwx(x,pso_start_x_fx_dir))
                            fx2 = fx
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write((' fx=%s\n' % fx))
                            xxxx.close()
                            breakx = xf(x, ('break_%s_%s-%s-%s_%s' % (numpam, id1, id2, id3, K)))
                            break

                    elif fx3 < fx2:
                        x1 = x2
                        fx1 = fx2
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write('\n' + (' x1=%s' % x1))
                        xxxx.close()
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write((' fx1=%s\n' % fx1))
                        xxxx.close()
                        x2m = {}

                        # m0 = 0.0001
                        m0 = x3 - 0.0005
                        xm0 = xf(m0, 'x2')
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write('\n' + (' x2=%s' % m0))
                        xxxx.close()
                        fm0 = float('%.6f'% f_pwx(m0,pso_start_x_fx_dir))
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write((' fx2=%s\n' % fm0))
                        xxxx.close()
                        x2m[m0] = fm0

                        # m1=0.1
                        m1 = float('%.4f'% (x1 + 0.1 * (x3 - x1)))
                        xm1 = xf(m1, 'x2')
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write('\n' + (' x2=%s' % m1))
                        xxxx.close()
                        fm1 = float('%.6f'% f_pwx(m1,pso_start_x_fx_dir))
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write((' fx2=%s\n' % fm1))
                        xxxx.close()
                        x2m[m1] = fm1

                        # m2=0.2
                        m2 = float('%.4f'% (x1 + 0.3 * (x3 - x1)))
                        xm2 = xf(m2, 'x2')
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write('\n' + (' x2=%s' % m2))
                        xxxx.close()
                        fm2 = float('%.6f'% f_pwx(m2,pso_start_x_fx_dir))
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write((' fx2=%s\n' % fm2))
                        xxxx.close()
                        x2m[m2] = fm2

                        # m3=0.3
                        m3 = float('%.4f'% (x1 + 0.5 * (x3 - x1)))
                        xm3 = xf(m3, 'x2')
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write('\n' + (' x2=%s' % m3))
                        xxxx.close()
                        fm3 = float('%.6f'% f_pwx(m3,pso_start_x_fx_dir))
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write((' fx2=%s\n' % fm3))
                        xxxx.close()
                        x2m[m3] = fm3

                        # m4 = 0.4
                        m4 = float('%.4f'% (x1 + 0.7 * (x3 - x1)))
                        xm4 = xf(m4, 'x2')
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write('\n' + (' x2=%s' % m4))
                        xxxx.close()
                        fm4 = float('%.6f'% f_pwx(m4,pso_start_x_fx_dir))
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write((' fx2=%s\n' % fm4))
                        xxxx.close()
                        x2m[m4] = fm4

                        # m5 = 0.5
                        m5 = float('%.4f'% (x1 + 0.9 * (x3 - x1)))
                        xm5 = xf(m5, 'x2')
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write('\n' + (' x2=%s' % m5))
                        xxxx.close()
                        fm5 = float('%.6f'% f_pwx(m5,pso_start_x_fx_dir))
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write((' fx2=%s\n' % fm5))
                        xxxx.close()
                        x2m[m5] = fm5
                        pso_start_x_fx_dir.update(x2m)
                        x2mmin = sorted(x2m.items(), key=lambda item: item[1])[0][0]
                        x2fmmin = sorted(x2m.items(), key=lambda item: item[1])[0][1]
                        if x2fmmin < fx3:
                            x2 = x2mmin
                            fx2 = x2fmmin
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write('\n' + (' x2=%s' % x2))
                            xxxx.close()
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write((' fx2=%s\n' % fx2))
                            xxxx.close()
                        else:
                            x = ffpam(xstart, [id1, id2, id3])
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write('\n' + (' x=%s' % x))
                            xxxx.close()
                            xx = xf(x, 'x')
                            fx = float('%.6f' % f_pwx(x,pso_start_x_fx_dir))
                            fx2 = fx
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write((' fx=%s\n' % fx))
                            xxxx.close()
                            breakx = xf(x, ('break_%s_%s-%s-%s_%s' % (numpam, id1, id2, id3, K)))
                            break

                if math.isinf(fx1) or fx1 > 2.766947558369898e+50:
                    fx1 = 1.0e+50
                if math.isinf(fx2) or fx2 > 2.766947558369898e+50:
                    fx2 = 1.0e+50
                if math.isinf(fx3) or fx3 > 2.766947558369898e+50:
                    fx3 = 1.0e+50
                print(fx1,fx2,fx3,x1,x2,x3)
                A = [[pow(x1, 2), x1, 1], [pow(x2, 2), x2, 1], [pow(x3, 2), x3, 1]]
                B = [fx1, fx2, fx3]
                X = np.linalg.solve(A, B)
                a0, a1, _ = X
                x = float('%.4f'% (- a1 / (2 * a0)))
                ee = abs(x2-x)
                xx = xf(x,'x')
                xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                xxxx.write('\n'+ (' x=%s' % x))
                xxxx.close()
                fx = float('%.6f'% f_pwx(x,pso_start_x_fx_dir))
                xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                xxxx.write((' fx=%s\n'% fx))
                xxxx.close()
                x_list.append(x)
                pso_start_x_fx_dir[x] = fx
                #
                if fx2 > fx:
                    if x2 > x:
                        x3 = x2
                        xx3 = xf(x2, 'x3')
                        fx3 = fx2
                        x2 = x
                        xx2 = xf(x, 'x2')
                        fx2 = fx
                    else:
                        x1 = x2
                        xx1 = xf(x2, 'x1')
                        fx1 = fx2
                        x2 = x
                        xx2 = xf(x, 'x2')
                        fx2 = fx
                else:
                    if x2 > x:
                        x1 = x
                        xx1 = xf(x, 'x1')
                        fx1 = fx
                    else:
                        x3 = x
                        xx3 = xf(x, 'x3')
                        fx3 = fx
                X2.append(x2)
                Y.append(fx2)
                print(' ：%f' % (k, x2))
                if k > nn or (len(x_list) > 10 and x_list[-1] == x_list[-2] == x_list[-3]) or (len(x_list) > 10 and x_list[-1] == x_list[-3] and x_list[-2] == x_list[-4]):
                    break

            if fx <= fx2:
                if fx > FFFF * 1.2 and FFFF != 0:
                    y = FFFF
                    x = ffpam(xstart, [id1, id2, id3])
                    X2.append(x)
                    Y.append(y)
                    # final_fun(X2, Y)
                    xstart = xf(x, 'start')
                    XX[numpam][K] = X2
                    YY[numpam][K] = Y
                else:
                    if x == x_start or x == x_end:
                        x = ffpam(xstart, [id1, id2, id3])
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write('\n' + (' x=%s' % x))
                        xxxx.close()
                        fx = fs0
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write((' fx=%s\n' % fx))
                        xxxx.close()
                    y = fx
                    print(' x:', x)
                    X2.append(x)
                    Y.append(y)
                    #final_fun(X2, Y)
                    xstart = xf(x, 'start')
                    XX[numpam][K] = X2
                    YY[numpam][K] = Y
            else:
                if fx2 > FFFF * 1.2 and FFFF != 0:

                    y = FFFF
                    x2 = ffpam(xstart,[id1,id2,id3])
                    X2.append(x2)
                    Y.append(y)
                    # final_fun(X2, Y)
                    xstart = xf(x2, 'start')
                    XX[numpam][K] = X2
                    YY[numpam][K] = Y

                else:
                    if x2 == x_start or x2 == x_end:
                        x2 = ffpam(xstart, [id1, id2, id3])
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write('\n' + (' x=%s' % x2))
                        xxxx.close()
                        fx2 = fs0
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write((' fx=%s\n' % fx2))
                        xxxx.close()
                    y = fx2
                    print(" x2", x2)
                    X2.append(x2)
                    Y.append(y)
                    #final_fun(X2, Y)
                    xstart = xf(x2, 'start')
                    XX[numpam][K] = X2
                    YY[numpam][K] = Y
            xwrite = open(os.getcwd() + os.sep + outdir + os.sep + 'xkpambest', 'a')
            x_x = XX[numpam][K][-1]
            f_x = YY[numpam][K][-1]
            bound = ((x_start, x_end),)
            if float(f_x) < min(float(f_a),float(f_b)):
                pso_x_fx_list = sorted(pso_start_x_fx_dir.items(), key=lambda item: item[1])
                x2 = s0
                fx2 = fs0
                for x_fx in pso_x_fx_list:
                    x = x_fx[0]
                    fx = x_fx[1]
                    if x_start < x < x_end:
                        if fx < fx2:
                            x2 = x
                            fx2 = fx
                xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                xxxx.write('\n The pso result x=%s' % (x2))
                xxxx.close()
                xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                xxxx.write(' fx=%s' % (fx2))
                xxxx.close()

                xstart = xf(x2, 'start')
                XX[numpam][K].append(x2)
                YY[numpam][K].append(fx2)

                error_flag = "          "
            elif float(f_a) == float(f_x) == float(f_b):
                x2 = s0
                fx2 = fs0
                xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                xxxx.write('\n The pso result x=%s' % (x2))
                xxxx.close()
                xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                xxxx.write(' fx=%s' % (fx2))
                xxxx.close()
                xstart = xf(x2, 'start')
                XX[numpam][K].append(x2)
                YY[numpam][K].append(fx2)
                error_flag = "                                =="

            elif float(f_a) < float(f_x) < float(f_b):
                pso_s = pso_scan(f, m_ffx, x_pam_start, var_name, bound,pso_start_x_fx_dir)
                pso_start_x_fx_dir.update(pso_s[1])
                pso_x_fx_list = sorted(pso_start_x_fx_dir.items(), key=lambda item: item[1])
                x2 = s0
                fx2 = fs0

                for x_fx in pso_x_fx_list:
                    x = x_fx[0]
                    fx = x_fx[1]
                    if x_start < x < x_end:
                        if fx < fx2:
                            x2 = x
                            fx2 = fx
                xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                xxxx.write('\n The pso result x=%s' % (x2))
                xxxx.close()
                xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                xxxx.write(' fx=%s' % (fx2))
                xxxx.close()

                xstart = xf(x2, 'start')
                XX[numpam][K].append(x2)
                YY[numpam][K].append(fx2)
                error_flag = "  up  " + " PSO_result\n" +  str(numpam).ljust(5) + str(id1) + '-' + str(id2) + '-' + str(id3).ljust(2) + \
                         '  x=' + str('%.4f' % x2).ljust(10) +\
                         '  y=' + str(fx2).ljust(15)
            elif float(f_a) > float(f_x) > float(f_b):
                pso_s = pso_scan(f, m_ffx, x_pam_start, var_name, bound,pso_start_x_fx_dir)
                pso_start_x_fx_dir.update(pso_s[1])
                pso_x_fx_list = sorted(pso_start_x_fx_dir.items(), key=lambda item: item[1])
                x2 = s0
                fx2 = fs0
                for x_fx in pso_x_fx_list:
                    x = x_fx[0]
                    fx = x_fx[1]
                    if x_start < x < x_end:
                        if fx < fx2:
                            x2 = x
                            fx2 = fx
                xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                xxxx.write('\n The pso result x=%s' % (x2))
                xxxx.close()
                xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                xxxx.write(' fx=%s' % (fx2))
                xxxx.close()

                xstart = xf(x2, 'start')
                XX[numpam][K].append(x2)
                YY[numpam][K].append(fx2)
                error_flag = "  down" + " PSO_result\n" +  str(numpam).ljust(5) + str(id1) + '-' + str(id2) + '-' + str(id3).ljust(2) + \
                         '  x=' + str('%.4f' % x2).ljust(10) +\
                         '  y=' + str(fx2).ljust(15)
            elif float(f_a) < float(f_x) and float(f_x) > float(f_b):
                pso_s = pso_scan(f, m_ffx, x_pam_start, var_name, bound,pso_start_x_fx_dir)
                pso_start_x_fx_dir.update(pso_s[1])
                pso_x_fx_list = sorted(pso_start_x_fx_dir.items(), key=lambda item: item[1])
                x2 = s0
                fx2 = fs0
                for x_fx in pso_x_fx_list:
                    x = x_fx[0]
                    fx = x_fx[1]
                    if x_start < x < x_end:
                        if fx < fx2:
                            x2 = x
                            fx2 = fx
                xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                xxxx.write('\n The pso result x=%s' % (x2))
                xxxx.close()
                xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                xxxx.write(' fx=%s' % (fx2))
                xxxx.close()

                xstart = xf(x2, 'start')
                XX[numpam][K].append(x2)
                YY[numpam][K].append(fx2)
                error_flag = "  凸" + " PSO_result\n" +  str(numpam).ljust(5) + str(id1) + '-' + str(id2) + '-' + str(id3).ljust(2) + \
                         '  x=' + str('%.4f' % x2).ljust(10) +\
                         '  y=' + str(fx2).ljust(15)
            else:
                pso_s = pso_scan(f, m_ffx, x_pam_start, var_name, bound,pso_start_x_fx_dir)
                pso_start_x_fx_dir.update(pso_s[1])
                pso_x_fx_list = sorted(pso_start_x_fx_dir.items(), key=lambda item: item[1])
                x2 = s0
                fx2 = fs0
                for x_fx in pso_x_fx_list:
                    x = x_fx[0]
                    fx = x_fx[1]
                    if x_start < x < x_end:
                        if fx < fx2:
                            x2 = x
                            fx2 = fx
                xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                xxxx.write('\n The pso result x=%s' % (x2))
                xxxx.close()
                xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                xxxx.write(' fx=%s' % (fx2))
                xxxx.close()

                xstart = xf(x2, 'start')
                XX[numpam][K].append(x2)
                YY[numpam][K].append(fx2)
                error_flag = "  not fit" + " PSO_result\n" +  str(numpam).ljust(5) + str(id1) + '-' + str(id2) + '-' + str(id3).ljust(2) + \
                         '  x=' + str('%.4f' % x2).ljust(10) +\
                         '  y=' + str(fx2).ljust(15)
            xwrite = open(os.getcwd() + os.sep + outdir + os.sep + 'xkpambest', 'a')
            xwrite.write(str(numpam).ljust(5) + str(id1) + '-' + str(id2) + '-' + str(id3).ljust(2) + \
                         '  x=' + str('%.4f' % x_x).ljust(10) +\
                         '  y=' + str(f_x).ljust(15) + \
                         ' a=' + str(x_a).ljust(10) + ' b=' + str(x_b).ljust(10) + \
                         ' fa=' + str(f_a).ljust(15) + ' fb=' + str(f_b).ljust(15) + error_flag + '\n')
            xwrite.close()
            FFFF = YY[numpam][K][-1]
            ###
            if abs(float(XX[numpam][K][-1]) - s0):
                re_num += 1
                re_pam[numpam] = {}
                re_pam[numpam]['id'] = {}
                re_pam[numpam]['value'] = {}
                re_pam[numpam]['id'][1] = id1
                re_pam[numpam]['id'][2] = id2
                re_pam[numpam]['id'][3] = id3
                re_pam[numpam]['value'][1] = a
                re_pam[numpam]['value'][2] = b
                re_pam[numpam]['step'] = pam[numpam]['step']
                with open(os.getcwd() + os.sep + outdir + os.sep + 're_params', 'a') as re:
                    re.write(str(id1).ljust(3) + str(id2).ljust(4) + str(id3).ljust(4) + "  0.1  "+\
                             str(a).rjust(12) + str(b).rjust(12) + '\n')
        bestwrite = ff(xstart, ('best%s' % K))
        if scale_flag == 0:
            write_params = w_params(dir_path + os.sep + outdir + os.sep + "ffield.reax.start")
            paramename = 'params'
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
        if re_num > 1:
            numpam = 0
            tt = {}
            ss = {}
            pam_final.update(re_pam)
            for key, val in re_pam.items():
                numpam = key
                k = 0
                a = pam[numpam]['value'][1]
                b = pam[numpam]['value'][2]
                X2, Y = list(), list()
                XX[numpam][K] = []
                YY[numpam][K] = []
                id1 = int(pam[numpam]['id'][1])
                id2 = int(pam[numpam]['id'][2])
                id3 = int(pam[numpam]['id'][3])
                pam_flag = str(pam[numpam]['id'][1]) + '-' + str(pam[numpam]['id'][2]) + '-' + str(pam[numpam]['id'][3])
                x2s = {}
                pso_start_x_fx_dir = {}
                #
                def xf(x, name):
                    global xstart
                    return ff(xstart, [id1, id2, id3, x, str(name)])

                #
                def f_pwx(x, pso_start_x_fx_dir):
                    pwx_x_name = xf(x, "pwx")
                    if x in pso_start_x_fx_dir:
                        f_v = pso_start_x_fx_dir[x]
                        print(" ")
                    else:
                        f_v = f(pwx_x_name)[0]
                        pso_start_x_fx_dir[x] = f_v
                    return f_v

                x_start = pam_scale[pam_flag]['a']
                x_end = pam_scale[pam_flag]['b']
                #
                x1 = float('%.4f' % a)
                x3 = float('%.4f' % b)
                # x2 = float('%.4f'% ((a + b) / 2.0))
                ee = e + 1
                xx1 = xf(x1, 'x1')
                # xx2 = xf(x2, 'x2')
                xx3 = xf(x3, 'x3')
                x_left = []
                x_right = []
                # x1
                xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                xxxx.write('\n' + str(K).ljust(3) + str(numpam).ljust(5) + str(id1).ljust(3) + str(id2).ljust(3) + str(
                    id3).ljust(3) + ' x1=%s' % (x1))
                xxxx.close()
                fx1 = float('%.6f' % f_pwx(x1, pso_start_x_fx_dir))
                xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                xxxx.write(' fx1=%s' % (fx1))
                xxxx.close()
                x_a = x1
                f_a = fx1
                x2s[x1] = fx1

                # x3
                xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                xxxx.write(' x3=%s' % (x3))
                xxxx.close()
                fx3 = float('%.6f' % f_pwx(x3, pso_start_x_fx_dir))
                xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                xxxx.write(' fx3=%s' % (fx3))
                xxxx.close()
                x_b = x3
                f_b = fx3
                x2s[x3] = fx3

                #
                s0 = ffpam(xstart, [id1, id2, id3])
                xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                xxxx.write('\n' + ' x0=%s' % (s0))
                xxxx.close()
                xs0 = xf(s0, 'x2')
                fs0 = float('%.6f' % f_pwx(s0, pso_start_x_fx_dir))
                xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                xxxx.write(' fx0=%s' % (fs0))
                xxxx.close()
                x2s[s0] = fs0
                x_left.append((x1, fx1))
                x_right.append((x3, fx3))

                x_pam_start = [s0]
                var_name = [pam_flag]
                bound = ((x_a, x_b),)
                #  pso
                pso_start_x_fx_dir.update(x2s)
                if mul_pso == 1:
                    pso_s = pso_scan(f, m_ffx, x_pam_start, var_name, bound, pso_start_x_fx_dir, num_pso=1, POP_SIZE=8)
                    pso_start_x_fx_dir.update(pso_s[1])

                    x2 = pso_s[0][pam_flag]
                    fx2 = pso_s[-1]
                    pso_x_fx_list = sorted(pso_s[1].items(), key=lambda item: item[1])
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write('\n The start x2=%s' % (x2))
                    xxxx.close()
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write(' fx2=%s                               PSO' % (fx2))
                    xxxx.close()
                    if x2 == x_start or x2 == x_end:
                        x2 = pso_x_fx_list[1][0]
                        fx2 = pso_x_fx_list[1][1]
                        pso_x_fx_list = sorted(pso_s[1].items(), key=lambda item: item[1])
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write('\n The start x2=%s' % (x2))
                        xxxx.close()
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write(' fx2=%s                               PSO' % (fx2))
                        xxxx.close()

                else:
                    if fs0 > fx1 or fs0 > fx3:
                        pso_s = pso_scan(f, m_ffx, x_pam_start, var_name, bound, pso_start_x_fx_dir, num_pso=1,
                                         POP_SIZE=8)
                        pso_start_x_fx_dir.update(pso_s[1])

                        x2 = pso_s[0][pam_flag]
                        fx2 = pso_s[-1]
                        pso_x_fx_list = sorted(pso_s[1].items(), key=lambda item: item[1])
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write('\n The start x2=%s' % (x2))
                        xxxx.close()
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write(' fx2=%s                               PSO' % (fx2))
                        xxxx.close()
                        if x2 == x_start or x2 == x_end:
                            x2 = pso_x_fx_list[1][0]
                            fx2 = pso_x_fx_list[1][1]
                            pso_x_fx_list = sorted(pso_s[1].items(), key=lambda item: item[1])
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write('\n The start x2=%s' % (x2))
                            xxxx.close()
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write(' fx2=%s                               PSO' % (fx2))
                            xxxx.close()
                    else:
                        x2 = s0
                        fx2 = fs0
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write('\n The start x2=%s' % (x2))
                        xxxx.close()
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write(' fx2=%s  ' % (fx2))
                        xxxx.close()
                if x2 == x_start or x2 == x_end:
                    x2 = s0
                    fx2 = fs0

                fx = fx2
                x = x2
                x_list = []
                x2m = {}
                while ee > e:
                    k += 1
                    #
                    if fx2 > fx1 or fx2 > fx3:
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + '凸点', 'a')
                        xxxx.write('\n' + ('%s %s %s %s %s %s %s %s %s' % (K, numpam, id1, id2, id3, k, fx1, fx2, fx3)))
                        xxxx.close()
                        #
                        fmmin = min(fx1, fx3)

                        # m0 = 0.0001
                        m0 = x1 + 0.0005
                        xm0 = xf(m0, 'x2')
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write('\n' + (' x2=%s' % m0))
                        xxxx.close()
                        fm0 = float('%.6f' % f_pwx(m0, pso_start_x_fx_dir))
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write((' fx2=%s\n' % fm0))
                        xxxx.close()
                        x2m[m0] = fm0

                        # m1 = 0.0125
                        m1 = float('%.4f' % (x1 + 0.0125 * (x3 - x1)))
                        xm1 = xf(m1, 'x2')
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write('\n' + (' x2=%s' % m1))
                        xxxx.close()
                        fm1 = float('%.6f' % f_pwx(m1, pso_start_x_fx_dir))
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write((' fx2=%s\n' % fm1))
                        xxxx.close()
                        x2m[m1] = fm1

                        # m3 = 0.05
                        m3 = float('%.4f' % (x1 + 0.05 * (x3 - x1)))
                        xm3 = xf(m3, 'x2')
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write('\n' + (' x2=%s' % m3))
                        xxxx.close()
                        fm3 = float('%.6f' % f_pwx(m3, pso_start_x_fx_dir))
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write((' fx2=%s\n' % fm3))
                        xxxx.close()
                        x2m[m3] = fm3

                        # m5 = 0.1
                        m5 = float('%.4f' % (x1 + 0.1 * (x3 - x1)))
                        xm5 = xf(m5, 'x2')
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write('\n' + (' x2=%s' % m5))
                        xxxx.close()
                        fm5 = float('%.6f' % f_pwx(m5, pso_start_x_fx_dir))
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write((' fx2=%s\n' % fm5))
                        xxxx.close()
                        x2m[m5] = fm5

                        # m7 = 0.9
                        m7 = float('%.4f' % (x1 + 0.9 * (x3 - x1)))
                        xm7 = xf(m7, 'x2')
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write('\n' + (' x2=%s' % m7))
                        xxxx.close()
                        fm7 = float('%.6f' % f_pwx(m7, pso_start_x_fx_dir))
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write((' fx2=%s\n' % fm7))
                        xxxx.close()
                        x2m[m7] = fm7

                        # m9 = 0.995
                        m9 = float('%.4f' % (x1 + 0.995 * (x3 - x1)))
                        xm9 = xf(m9, 'x2')
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write('\n' + (' x2=%s' % m9))
                        xxxx.close()
                        fm9 = float('%.6f' % f_pwx(m9, pso_start_x_fx_dir))
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                        xxxx.write((' fx2=%s\n' % fm9))
                        xxxx.close()
                        x2m[m9] = fm9
                        pso_start_x_fx_dir.update(x2m)
                        x2mmin = sorted(x2m.items(), key=lambda item: item[1])[0][0]
                        x2fmmin = sorted(x2m.items(), key=lambda item: item[1])[0][1]
                        if x2fmmin < fmmin:
                            x2 = x2mmin
                            fx2 = x2fmmin
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write('\n' + (' x2=%s' % x2))
                            xxxx.close()
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write((' fx2=%s\n' % fx2))
                            xxxx.close()
                        else:

                            x = ffpam(xstart, [id1, id2, id3])
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write('\n' + (' x=%s' % x))
                            xxxx.close()
                            xx = xf(x, 'x')
                            fx = float('%.6f' % f_pwx(x, pso_start_x_fx_dir))
                            fx2 = fx
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write((' fx=%s\n' % fx))
                            xxxx.close()
                            breakx = xf(x, ('break_%s_%s-%s-%s_%s' % (numpam, id1, id2, id3, K)))
                            break

                    #
                    elif k == 1 and (fx1 == fx2 or fx2 == fx3):
                        xxxx = open(os.getcwd() + os.sep + outdir + os.sep + ' ', 'a')
                        xxxx.write('\n' + ('%s %s %s %s %s %s %s %s %s' % (K, numpam, id1, id2, id3, k, fx1, fx2, fx3)))
                        xxxx.close()
                        if fx1 == fx2 == fx3:
                            x = ffpam(xstart, [id1, id2, id3])
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write('\n' + (' x=%s' % x))
                            xxxx.close()
                            fx = fx1
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write((' fx=%s\n' % fx))
                            xxxx.close()
                            print(" ")
                            breakx3 = xf(x1, ('break_%s_%s-%s-%s_%s' % (numpam, id1, id2, id3, K)))
                            break

                        elif fx1 > fx2:
                            x = x2
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write('\n' + (' x=%s' % x))
                            xxxx.close()

                            fx = fx2
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write((' fx=%s\n' % fx))
                            xxxx.close()
                            breakx = xf(x, ('break_%s_%s-%s-%s_%s' % (numpam, id1, id2, id3, K)))
                            break
                        elif fx3 > fx2:
                            x = x2
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write('\n' + (' x=%s' % x))
                            xxxx.close()

                            fx = fx2
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write((' fx=%s\n' % fx))
                            xxxx.close()
                            breakx = xf(x, ('break_%s_%s-%s-%s_%s' % (numpam, id1, id2, id3, K)))
                            break
                        elif fx1 < fx2:
                            x3 = x2
                            fx3 = fx2
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write('\n' + (' x3=%s' % x3))
                            xxxx.close()
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write((' fx3=%s\n' % fx3))
                            xxxx.close()
                            #
                            x2m = {}

                            # m0 = 0.0001
                            m0 = x1 + 0.0005
                            xm0 = xf(m0, 'x2')
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write('\n' + (' x2=%s' % m0))
                            xxxx.close()
                            fm0 = float('%.6f' % f_pwx(m0, pso_start_x_fx_dir))
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write((' fx2=%s\n' % fm0))
                            xxxx.close()
                            x2m[m0] = fm0

                            # m1=0.1
                            m1 = float('%.4f' % (x1 + 0.1 * (x3 - x1)))
                            xm1 = xf(m1, 'x2')
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write('\n' + (' x2=%s' % m1))
                            xxxx.close()
                            fm1 = float('%.6f' % f_pwx(m1, pso_start_x_fx_dir))
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write((' fx2=%s\n' % fm1))
                            xxxx.close()
                            x2m[m1] = fm1

                            # m2=0.2
                            m2 = float('%.4f' % (x1 + 0.2 * (x3 - x1)))
                            xm2 = xf(m2, 'x2')
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write('\n' + (' x2=%s' % m2))
                            xxxx.close()
                            fm2 = float('%.6f' % f_pwx(m2, pso_start_x_fx_dir))
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write((' fx2=%s\n' % fm2))
                            xxxx.close()
                            x2m[m2] = fm2

                            # m3=0.3
                            m3 = float('%.4f' % (x1 + 0.3 * (x3 - x1)))
                            xm3 = xf(m3, 'x2')
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write('\n' + (' x2=%s' % m3))
                            xxxx.close()
                            fm3 = float('%.6f' % f_pwx(m3, pso_start_x_fx_dir))
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write((' fx2=%s\n' % fm3))
                            xxxx.close()
                            x2m[m3] = fm3

                            # m4 = 0.4
                            m4 = float('%.4f' % (x1 + 0.4 * (x3 - x1)))
                            xm4 = xf(m4, 'x2')
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write('\n' + (' x2=%s' % m4))
                            xxxx.close()
                            fm4 = float('%.6f' % f_pwx(m4, pso_start_x_fx_dir))
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write((' fx2=%s\n' % fm4))
                            xxxx.close()
                            x2m[m4] = fm4

                            # m5 = 0.5
                            m5 = float('%.4f' % (x1 + 0.5 * (x3 - x1)))
                            xm5 = xf(m5, 'x2')
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write('\n' + (' x2=%s' % m5))
                            xxxx.close()
                            fm5 = float('%.6f' % f_pwx(m5, pso_start_x_fx_dir))
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write((' fx2=%s\n' % fm5))
                            xxxx.close()
                            x2m[m5] = fm5
                            pso_start_x_fx_dir.update(x2m)
                            x2mmin = sorted(x2m.items(), key=lambda item: item[1])[0][0]
                            x2fmmin = sorted(x2m.items(), key=lambda item: item[1])[0][1]
                            if x2fmmin < fx1:
                                x2 = x2mmin
                                fx2 = x2fmmin
                                xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                                xxxx.write('\n' + (' x2=%s' % x2))
                                xxxx.close()
                                xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                                xxxx.write((' fx2=%s\n' % fx2))
                                xxxx.close()
                            else:
                                x = ffpam(xstart, [id1, id2, id3])
                                xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                                xxxx.write('\n' + (' x=%s' % x))
                                xxxx.close()
                                xx = xf(x, 'x')
                                fx = float('%.6f' % f_pwx(x, pso_start_x_fx_dir))
                                fx2 = fx
                                xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                                xxxx.write((' fx=%s\n' % fx))
                                xxxx.close()
                                breakx = xf(x, ('break_%s_%s-%s-%s_%s' % (numpam, id1, id2, id3, K)))
                                break

                        elif fx3 < fx2:
                            x1 = x2
                            fx1 = fx2
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write('\n' + (' x1=%s' % x1))
                            xxxx.close()
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write((' fx1=%s\n' % fx1))
                            xxxx.close()
                            #
                            x2m = {}

                            # m0 = 0.0001
                            m0 = x3 - 0.0005
                            xm0 = xf(m0, 'x2')
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write('\n' + (' x2=%s' % m0))
                            xxxx.close()
                            fm0 = float('%.6f' % f_pwx(m0, pso_start_x_fx_dir))
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write((' fx2=%s\n' % fm0))
                            xxxx.close()
                            x2m[m0] = fm0

                            # m1=0.1
                            m1 = float('%.4f' % (x1 + 0.1 * (x3 - x1)))
                            xm1 = xf(m1, 'x2')
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write('\n' + (' x2=%s' % m1))
                            xxxx.close()
                            fm1 = float('%.6f' % f_pwx(m1, pso_start_x_fx_dir))
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write((' fx2=%s\n' % fm1))
                            xxxx.close()
                            x2m[m1] = fm1

                            # m2=0.2
                            m2 = float('%.4f' % (x1 + 0.3 * (x3 - x1)))
                            xm2 = xf(m2, 'x2')
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write('\n' + (' x2=%s' % m2))
                            xxxx.close()
                            fm2 = float('%.6f' % f_pwx(m2, pso_start_x_fx_dir))
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write((' fx2=%s\n' % fm2))
                            xxxx.close()
                            x2m[m2] = fm2

                            # m3=0.3
                            m3 = float('%.4f' % (x1 + 0.5 * (x3 - x1)))
                            xm3 = xf(m3, 'x2')
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write('\n' + (' x2=%s' % m3))
                            xxxx.close()
                            fm3 = float('%.6f' % f_pwx(m3, pso_start_x_fx_dir))
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write((' fx2=%s\n' % fm3))
                            xxxx.close()
                            x2m[m3] = fm3

                            # m4 = 0.4
                            m4 = float('%.4f' % (x1 + 0.7 * (x3 - x1)))
                            xm4 = xf(m4, 'x2')
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write('\n' + (' x2=%s' % m4))
                            xxxx.close()
                            fm4 = float('%.6f' % f_pwx(m4, pso_start_x_fx_dir))
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write((' fx2=%s\n' % fm4))
                            xxxx.close()
                            x2m[m4] = fm4

                            # m5 = 0.5
                            m5 = float('%.4f' % (x1 + 0.9 * (x3 - x1)))
                            xm5 = xf(m5, 'x2')
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write('\n' + (' x2=%s' % m5))
                            xxxx.close()
                            fm5 = float('%.6f' % f_pwx(m5, pso_start_x_fx_dir))
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write((' fx2=%s\n' % fm5))
                            xxxx.close()
                            x2m[m5] = fm5
                            pso_start_x_fx_dir.update(x2m)
                            x2mmin = sorted(x2m.items(), key=lambda item: item[1])[0][0]
                            x2fmmin = sorted(x2m.items(), key=lambda item: item[1])[0][1]
                            if x2fmmin < fx3:
                                x2 = x2mmin
                                fx2 = x2fmmin
                                xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                                xxxx.write('\n' + (' x2=%s' % x2))
                                xxxx.close()
                                xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                                xxxx.write((' fx2=%s\n' % fx2))
                                xxxx.close()
                            else:
                                x = ffpam(xstart, [id1, id2, id3])
                                xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                                xxxx.write('\n' + (' x=%s' % x))
                                xxxx.close()
                                xx = xf(x, 'x')
                                fx = float('%.6f' % f_pwx(x, pso_start_x_fx_dir))
                                fx2 = fx
                                xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                                xxxx.write((' fx=%s\n' % fx))
                                xxxx.close()
                                breakx = xf(x, ('break_%s_%s-%s-%s_%s' % (numpam, id1, id2, id3, K)))
                                break

                    if math.isinf(fx1) or fx1 > 2.766947558369898e+50:
                        fx1 = 1.0e+50
                    if math.isinf(fx2) or fx2 > 2.766947558369898e+50:
                        fx2 = 1.0e+50
                    if math.isinf(fx3) or fx3 > 2.766947558369898e+50:
                        fx3 = 1.0e+50
                    print(fx1, fx2, fx3, x1, x2, x3)
                    A = [[pow(x1, 2), x1, 1], [pow(x2, 2), x2, 1], [pow(x3, 2), x3, 1]]
                    B = [fx1, fx2, fx3]
                    X = np.linalg.solve(A, B)
                    a0, a1, _ = X
                    x = float('%.4f' % (- a1 / (2 * a0)))
                    ee = abs(x2 - x)
                    xx = xf(x, 'x')
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write('\n' + (' x=%s' % x))
                    xxxx.close()
                    fx = float('%.6f' % f_pwx(x, pso_start_x_fx_dir))
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write((' fx=%s\n' % fx))
                    xxxx.close()
                    x_list.append(x)
                    pso_start_x_fx_dir[x] = fx

                    if fx2 > fx:
                        if x2 > x:
                            x3 = x2
                            xx3 = xf(x2, 'x3')
                            fx3 = fx2
                            x2 = x
                            xx2 = xf(x, 'x2')
                            fx2 = fx
                        else:
                            x1 = x2
                            xx1 = xf(x2, 'x1')
                            fx1 = fx2
                            x2 = x
                            xx2 = xf(x, 'x2')
                            fx2 = fx
                    else:
                        if x2 > x:
                            x1 = x
                            xx1 = xf(x, 'x1')
                            fx1 = fx
                        else:
                            x3 = x
                            xx3 = xf(x, 'x3')
                            fx3 = fx
                    X2.append(x2)
                    Y.append(fx2)
                    print('  ：%f' % (k, x2))
                    if k > nn or (len(x_list) > 10 and x_list[-1] == x_list[-2] == x_list[-3]) or (
                            len(x_list) > 10 and x_list[-1] == x_list[-3] and x_list[-2] == x_list[-4]):
                        break
                    # update_point(X2, Y)

                if fx <= fx2:
                    if fx > FFFF * 1.2 and FFFF != 0:
                        y = FFFF
                        x = ffpam(xstart, [id1, id2, id3])
                        X2.append(x)
                        Y.append(y)
                        # final_fun(X2, Y)
                        xstart = xf(x, 'start')
                        XX[numpam][K] = X2
                        YY[numpam][K] = Y
                    else:
                        if x == x_start or x == x_end:
                            x = ffpam(xstart, [id1, id2, id3])
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write('\n' + (' x=%s' % x))
                            xxxx.close()
                            fx = fs0
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write((' fx=%s\n' % fx))
                            xxxx.close()
                        y = fx
                        print(' x:', x)
                        X2.append(x)
                        Y.append(y)
                        # final_fun(X2, Y)
                        xstart = xf(x, 'start')
                        XX[numpam][K] = X2
                        YY[numpam][K] = Y
                else:
                    if fx2 > FFFF * 1.2 and FFFF != 0:

                        y = FFFF
                        x2 = ffpam(xstart, [id1, id2, id3])
                        X2.append(x2)
                        Y.append(y)
                        # final_fun(X2, Y)
                        xstart = xf(x2, 'start')
                        XX[numpam][K] = X2
                        YY[numpam][K] = Y

                    else:
                        if x2 == x_start or x2 == x_end:
                            x2 = ffpam(xstart, [id1, id2, id3])
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write('\n' + (' x=%s' % x2))
                            xxxx.close()
                            fx2 = fs0
                            xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                            xxxx.write((' fx=%s\n' % fx2))
                            xxxx.close()
                        y = fx2
                        print(" x2", x2)
                        X2.append(x2)
                        Y.append(y)
                        # final_fun(X2, Y)
                        xstart = xf(x2, 'start')
                        XX[numpam][K] = X2
                        YY[numpam][K] = Y
                xwrite = open(os.getcwd() + os.sep + outdir + os.sep + 'xkpambest', 'a')
                x_x = XX[numpam][K][-1]
                f_x = YY[numpam][K][-1]
                ############################  PSO scan
                bound = ((x_start, x_end),)
                if float(f_x) < min(float(f_a), float(f_b)):
                    pso_x_fx_list = sorted(pso_start_x_fx_dir.items(), key=lambda item: item[1])
                    x2 = s0
                    fx2 = fs0
                    for x_fx in pso_x_fx_list:
                        x = x_fx[0]
                        fx = x_fx[1]
                        if x_start < x < x_end:
                            if fx < fx2:
                                x2 = x
                                fx2 = fx
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write('\n The pso result x=%s' % (x2))
                    xxxx.close()
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write(' fx=%s' % (fx2))
                    xxxx.close()

                    xstart = xf(x2, 'start')
                    XX[numpam][K].append(x2)
                    YY[numpam][K].append(fx2)

                    error_flag = "          "
                elif float(f_a) == float(f_x) == float(f_b):
                    x2 = s0
                    fx2 = fs0
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write('\n The pso result x=%s' % (x2))
                    xxxx.close()
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write(' fx=%s' % (fx2))
                    xxxx.close()
                    xstart = xf(x2, 'start')
                    XX[numpam][K].append(x2)
                    YY[numpam][K].append(fx2)
                    error_flag = "                                =="
                elif float(f_a) < float(f_x) < float(f_b):
                    pso_s = pso_scan(f, m_ffx, x_pam_start, var_name, bound, pso_start_x_fx_dir)
                    pso_start_x_fx_dir.update(pso_s[1])
                    pso_x_fx_list = sorted(pso_start_x_fx_dir.items(), key=lambda item: item[1])
                    x2 = s0
                    fx2 = fs0

                    for x_fx in pso_x_fx_list:
                        x = x_fx[0]
                        fx = x_fx[1]
                        if x_start < x < x_end:
                            if fx < fx2:
                                x2 = x
                                fx2 = fx
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write('\n The pso result x=%s' % (x2))
                    xxxx.close()
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write(' fx=%s' % (fx2))
                    xxxx.close()

                    xstart = xf(x2, 'start')
                    XX[numpam][K].append(x2)
                    YY[numpam][K].append(fx2)
                    error_flag = "  up  " + " PSO_result\n" + str(numpam).ljust(5) + str(id1) + '-' + str(
                        id2) + '-' + str(id3).ljust(2) + \
                                 '  x=' + str('%.4f' % x2).ljust(10) + \
                                 '  y=' + str(fx2).ljust(15)
                elif float(f_a) > float(f_x) > float(f_b):
                    pso_s = pso_scan(f, m_ffx, x_pam_start, var_name, bound, pso_start_x_fx_dir)
                    pso_start_x_fx_dir.update(pso_s[1])
                    pso_x_fx_list = sorted(pso_start_x_fx_dir.items(), key=lambda item: item[1])
                    x2 = s0
                    fx2 = fs0
                    for x_fx in pso_x_fx_list:
                        x = x_fx[0]
                        fx = x_fx[1]
                        if x_start < x < x_end:
                            if fx < fx2:
                                x2 = x
                                fx2 = fx
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write('\n The pso result x=%s' % (x2))
                    xxxx.close()
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write(' fx=%s' % (fx2))
                    xxxx.close()

                    xstart = xf(x2, 'start')
                    XX[numpam][K].append(x2)
                    YY[numpam][K].append(fx2)
                    error_flag = "  down" + " PSO_result\n" + str(numpam).ljust(5) + str(id1) + '-' + str(
                        id2) + '-' + str(id3).ljust(2) + \
                                 '  x=' + str('%.4f' % x2).ljust(10) + \
                                 '  y=' + str(fx2).ljust(15)
                elif float(f_a) < float(f_x) and float(f_x) > float(f_b):
                    pso_s = pso_scan(f, m_ffx, x_pam_start, var_name, bound, pso_start_x_fx_dir)
                    pso_start_x_fx_dir.update(pso_s[1])
                    pso_x_fx_list = sorted(pso_start_x_fx_dir.items(), key=lambda item: item[1])
                    x2 = s0
                    fx2 = fs0
                    for x_fx in pso_x_fx_list:
                        x = x_fx[0]
                        fx = x_fx[1]
                        if x_start < x < x_end:
                            if fx < fx2:
                                x2 = x
                                fx2 = fx
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write('\n The pso result x=%s' % (x2))
                    xxxx.close()
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write(' fx=%s' % (fx2))
                    xxxx.close()

                    xstart = xf(x2, 'start')
                    XX[numpam][K].append(x2)
                    YY[numpam][K].append(fx2)
                    error_flag = "   " + " PSO_result\n" + str(numpam).ljust(5) + str(id1) + '-' + str(id2) + '-' + str(
                        id3).ljust(2) + \
                                 '  x=' + str('%.4f' % x2).ljust(10) + \
                                 '  y=' + str(fx2).ljust(15)
                else:
                    pso_s = pso_scan(f, m_ffx, x_pam_start, var_name, bound, pso_start_x_fx_dir)
                    pso_start_x_fx_dir.update(pso_s[1])
                    pso_x_fx_list = sorted(pso_start_x_fx_dir.items(), key=lambda item: item[1])
                    x2 = s0
                    fx2 = fs0
                    for x_fx in pso_x_fx_list:
                        x = x_fx[0]
                        fx = x_fx[1]
                        if x_start < x < x_end:
                            if fx < fx2:
                                x2 = x
                                fx2 = fx
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write('\n The pso result x=%s' % (x2))
                    xxxx.close()
                    xxxx = open(os.getcwd() + os.sep + outdir + os.sep + 'xxxx', 'a')
                    xxxx.write(' fx=%s' % (fx2))
                    xxxx.close()

                    xstart = xf(x2, 'start')
                    XX[numpam][K].append(x2)
                    YY[numpam][K].append(fx2)
                    error_flag = "  not fit" + " PSO_result\n" + str(numpam).ljust(5) + str(id1) + '-' + str(
                        id2) + '-' + str(id3).ljust(2) + \
                                 '  x=' + str('%.4f' % x2).ljust(10) + \
                                 '  y=' + str(fx2).ljust(15)
                xwrite = open(os.getcwd() + os.sep + outdir + os.sep + 'xkpambest', 'a')
                xwrite.write(str(numpam).ljust(5) + str(id1) + '-' + str(id2) + '-' + str(id3).ljust(2) + \
                             '  x=' + str('%.4f' % x_x).ljust(10) + \
                             '  y=' + str(f_x).ljust(15) + \
                             ' a=' + str(x_a).ljust(10) + ' b=' + str(x_b).ljust(10) + \
                             ' fa=' + str(f_a).ljust(15) + ' fb=' + str(f_b).ljust(15) + error_flag + '\n')
                xwrite.close()
                FFFF = YY[numpam][K][-1]
            bestwrite = ff(xstart, ('best%s' % (K + 100)))

        EF[K] = YY[numpam][K][-1]
        error = open(os.getcwd() + os.sep + outdir + os.sep + 'ferror', 'a')
        error.write(str(K).ljust(5) + str(EF[K]) + '\n')
        error.close()
        if (EF[K-1] - EF[K]) < 0:
            nK += 1
            NK[nK] = K
            dertaK = NK[nK] - NK[nK-1]
        #####################
        if scale_flag == 0:
            write_params = w_params(dir_path + os.sep + outdir + os.sep + "ffield.reax.start")
            paramename = 'params'
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
    wbest = ff(xstart, ['best'])
    print("fffff：%d " % (K - 1))

    return(xstart,K,XX,YY,EF)
