# -*- coding:utf-8 -*-
import sys
def get_argv(argv):
    a_ = "    \"-a\"  ((QM-ff)/δ) * ((QM-ff)/δ);(defauld)\n\n"
    b_ = "    \"-b\"  ((QM-ff)/δ) * ((QM-ff)/δ) * weight\n\
              weight = entry-weight (wi)\n\n"
    c_ = "    \"-c\"  ((QM-ff/δ) * (QM-ff/δ)) * weight * 100,\n\
              weight = (wi/sum-sec-wi) * sec-weight,\n\
              sec-weight = num-sec-entry /num-all-entry\n\n"
    d_ = "    \"-d\"  followed by 5 num: w1,w2,w3,w4,w5;\n\
              w1 charge, w2 geometry, w3 energy, w4=0 cell-parameter, w5=0 heat\n\
              ((QM-ff/δ) * (QM-ff/δ)) * weight * 100,\n\
              weight = (wi/sum-sec-wi) * sec-weight,\n\
              sec-weight = wi /(w1+w2+w3+w4)\n\n"
    n_ = "The sub-folder \"data\" used for lammps-data file\n\
    The \"outdata\" folder used for output file\n\nThe \"ffield.reax.best\" form files is best ffiled for just 1-big step \n\n\
    The \"ffield.reax.start\" form file is the best ffield for every steps\n     can be used as restart file named \"ffield\" "

    ##############################
    object_function = {"-a":a_,"-b":b_,"-c":c_,"-d":d_}
    function_flag   = "-a"
    mpi_nc,mnp,md,m ,scale_flag= 4,2,0,0,0
    single_scan_flag,mul_pso,pso_flag,pso_new_flag,sa_flag,ga_flag,paowu_plus_pso_flag,mini_all_flag = 1,0,0,0,0,0,0,0
    fc, fg, fe, fl, fh, fm = 2, 2, 2, 2, 2, 2
    ###############################             help
    if "-s" in argv:
        scale_flag = 1
    if "-np" in argv:
        try:
            mpi_nc = int(argv[argv.index("-np") + 1])
        except:
            print("Command-line error:-np must followed by  number")
            sys.exit()
    if "-mnp" in argv:
        try:
            mnp = int(argv[argv.index("-mnp") + 1])
        except:
            print("Command-line error:-mnp must followed by  number")
            sys.exit()
    if "-m" in argv:
        print("minimize all geo")
        mini_all_flag = 1
        m = 1
    if "-a" in argv or "-b" in argv or "-c" in argv or "-d" in argv:
        if "-a" in argv:
            fc = fg = fe = fl = fh = fm =2
            function_flag = "-a"
        elif "-b" in argv:
            fc = fg = fe = fl = fh = fm =3
            function_flag = "-b"
        elif "-c" in argv:
            fc = fg = fe = fl = fh = fm =1
            function_flag = "-c"
        elif "-d" in argv:
            if len(argv) < 7:
                print("Command-line error: \"-d\" must followed by 6 numbers")
                sys.exit()
            argv_d = []
            try:
                argv_d = argv[argv.index("-d")+1 : argv.index("-d") + 7]
                fc = float(argv_d[0])
                fg = float(argv_d[1])
                fe = float(argv_d[2])
                fl = float(argv_d[3])
                fh = float(argv_d[4])
                fm = float(argv_d[5])
            except:
                print("Command-line error: \"-d\" must be followed by 6 numbers")
                sys.exit()
            else:
                function_flag = "-d"
                print(fc,fg,fe,fl,fh,fm)
    print("\n\nUse " + str(mpi_nc) + " threads\n\n")
    print("Use object function form as :\n" + object_function[function_flag])
    print("    Use command-line to change the object function:\n")
    print(a_ + b_ + c_ + d_ + n_)
    #####  Fix error
    fix_error = {}
    if "-f" in argv:
        try:
            fix_file = argv[argv.index("-f") + 1]
            fix_sum = 0
            fix_num = []
            with open(os.getcwd() + os.sep + fix_file) as r_fix:
                fix_lines = r_fix.readlines()
                fix_flag = 0
                start_num = 0
                end_num = 0
                for line in fix_lines:
                    spline = line.split()
                    if spline[0].isdigit():
                        job_num = int(spline[0])
                        if spline[-1] == "fix_start":
                            fix_flag = 1
                            start_num += 1
                        elif spline[-1] == "fix_end":
                            fix_flag = 0
                            end_num += 1
                            if end_num != start_num:
                                sys.exit()
                        if fix_flag == 1:
                            fix_sum += 1
                            if spline[-1] == "fix_start":
                                fix_error[job_num] = float(spline[-4])
                            else:
                                fix_error[job_num] = float(spline[-3])
            with open(os.getcwd() + os.sep + fix_file) as r_fix:
                fix_lines = r_fix.readlines()
                for line in fix_lines:
                    spline = line.split()
                    if spline[0].isdigit():
                        job_num = int(spline[0])
                        if spline[-1] == "fix":
                            fix_sum += 1
                            fix_num.append(job_num)
                            fix_error[job_num] = float(spline[-4])
                        else:
                            pass
            if fix_sum == 0:
                print("There is no entryes to fix,please use \"fix\" as a flag to fix at the end of the line")
                sys.exit()
            else:
                print("\n\nfix entrys num. as:",fix_num)
        except:
            print("\nCommand-line error: \"-f\" must be followed by fix_error file")
            sys.exit()
    #####
    if "-md" in argv:
        md = 1
    if "-mul" in argv:
        mul_pso = 1
    if "-pso" in argv:
        pso_flag = 1
        single_scan_flag = 0
    elif "-ga" in argv:
        ga_flag = 1
        single_scan_flag = 0
        pso_flag = 0
    elif "-pso_new" in argv:
        pso_new_flag = 1
        ga_flag = 0
        single_scan_flag = 0
        pso_flag = 0
    elif "-sa" in argv:
        sa_flag = 1
        pso_new_flag = 0
        ga_flag = 0
        single_scan_flag = 0
        pso_flag = 0
    elif "-sp" in argv:
        pso_flag = 0
        pso_new_flag = 0
        ga_flag = 0
        single_scan_flag = 0
        paowu_plus_pso_flag = 1
    algorithm_flag = [single_scan_flag, mul_pso, pso_flag, pso_new_flag, sa_flag, ga_flag, paowu_plus_pso_flag,
                      mini_all_flag]
    obj_func_flag = [fc, fg, fe, fl, fh, fm]
    mpi_num = [mpi_nc, mnp, md, m,scale_flag]
    print(mpi_num ,algorithm_flag,obj_func_flag)
    return mpi_num ,algorithm_flag,obj_func_flag,fix_error
if __name__ == '__main__':
    argv = sys.argv
    print(argv)
    get_argv(argv)



