import os, sys
import string
from typing import List
argv = sys.argv
def w_params(ff_path = "ffield"):
    if "-scale" in argv:
        scale = float(argv[argv.index("-scale") + 1])
    else:
        scale = 1
    ffield_path = ff_path
    path = os.getcwd()
    paramename = path + os.sep + 'params.scale'
    openparam = open(paramename, 'r')  #
    pamlines = openparam.readlines()
    openparam.close()
    params_start = {}
    for i in pamlines:
        perpam = i.split()
        if len(perpam) > 0 and i.find("#") != 0:
            id1 = perpam[0]
            id2 = perpam[1]
            id3 = perpam[2]
            id_flag = id1 + '-' + id2 + '-' + id3
            params_start[id_flag] = {}
            params_start[id_flag]["a"] = float(perpam[-2])
            params_start[id_flag]["b"] = float(perpam[-1])
    ff = open(ffield_path, 'r')
    fflines = ff.readlines()
    ffield = {}

    i = 1
    j = 1
    k = 1
    value = 0.0
    numi = 0
    numj = 0
    numk = 0
    sectitle = []

    numsec1 =0
    numsec2 =0
    numsec3 =0
    numsec4 =0
    numsec5 =0
    numsec6 =0
    numsec7 =0

    sec1ids = 0
    sec2ids = 0
    sec3ids = 0
    sec4ids = 0
    sec5ids = 0
    sec6ids = 0
    sec7ids = 0

    sec1tid = 0
    sec2tid = 0
    sec3tid = 0
    sec4tid = 0
    sec5tid = 0
    sec6tid = 0
    sec7tid = 0

    atomid = []
    int2id = []
    int3id = []
    int4id = []
    lg = 0

    atomtypes: List[str] = ['C','H','O','N','Al','Si','S','Cl','X','Cu','Pt','Zr','Ni','Co','He','Ne','Ar','Kr','Xe','F']

    #
    for line in fflines:
        spline = line.split()
        if ('!' in line) and (spline[0].isdigit()):
            numi += 1
            ffield[numi] = {}
            sectitle.append(line)
        if spline[0] in atomtypes:
            bondt = fflines.index(line)
            if bondt >39:
                atomid.append(bondt)
        if spline[0].isdigit() and spline[1].isdigit():
            int2id.append(fflines.index(line))
        if spline[0].isdigit() and spline[1].isdigit() and spline[2].isdigit():
            int3id.append(fflines.index(line))
        if spline[0].isdigit() and spline[1].isdigit() and spline[2].isdigit() and spline[3].isdigit():
            int4id.append(fflines.index(line))
        else:
            pass

    numsec1 = int(sectitle[0].split()[0])
    numsec2 = int(sectitle[1].split()[0])
    numsec3 = int(sectitle[2].split()[0])
    numsec4 = int(sectitle[3].split()[0])
    numsec5 = int(sectitle[4].split()[0])
    numsec6 = int(sectitle[5].split()[0])
    numsec7 = int(sectitle[6].split()[0])

    sec1tid = fflines.index(sectitle[0])
    sec2tid = fflines.index(sectitle[1])
    sec3tid = fflines.index(sectitle[2])
    sec4tid = fflines.index(sectitle[3])
    sec5tid = fflines.index(sectitle[4])
    sec6tid = fflines.index(sectitle[5])
    sec7tid = fflines.index(sectitle[6])

    sec1ids = sec1tid + 1
    sec2ids = atomid[0]
    sec3ids = int2id[0]
    sec4ids = sec4tid + 1
    sec5ids = int3id[0]
    sec6ids = int4id[0]
    sec7ids = sec7tid + 1


    if atomid[1]-atomid[0] == 5 :
        lg = 1

    # sec 1

    numi = 0
    numj = 0
    numk = 0
    ffield[1][0] = {}
    ffi = {}
    for line in fflines:
        spline = line.split()
        if ((('!' in line) and type(eval(line.split()[0])) == float)):
            numj += 1
            ffield[1][0][numj] = spline[0]
            id_flag = "1-0-" + str(numj)
            ff_value = float(spline[0])
            ffi[id_flag] = ff_value

    # sec2

    numi = 0
    numj = 1
    numk = 0
    for ii in range(numsec2):
        ffield[2][(ii) + 1] = {}
        numk = 0
        for line in fflines[(sec2ids + ii * (4 + lg)): ((sec2ids + ii * (4 + lg)) + (4 + lg))]:
            for value in line.split():
                numk += 1
                ffield[2][(ii + 1)][numk] = line.split()[line.split().index(value)]
                if numk > 1:
                    id_flag = "2-" + str(ii + 1) + '-' + str(numk - 1)
                    ff_value = float(line.split()[line.split().index(value)])
                    ffi[id_flag] = ff_value

    # sec3

    numi = 0
    numj = 1
    numk = 0
    for ii in range(numsec3):
        ffield[3][(ii + 1)] = {}
        numk = 0
        for line in fflines[(sec3ids + ii * 2): ((sec3ids + ii * 2) + 2)]:
            for value in line.split():
                numk += 1
                ffield[3][((ii + 1))][numk] = line.split()[line.split().index(value)]
                if numk > 2:
                    id_flag = "3-" + str(ii + 1) + '-' + str(numk - 2)
                    ff_value = float(line.split()[line.split().index(value)])
                    ffi[id_flag] = ff_value

    # sec 4
    numi = 0
    numj = 1
    numk = 0
    for ii in range(numsec4):
        ffield[4][(ii + 1)] = {}
        numk = 0
        for line in fflines[(sec4ids + ii): (sec4ids + ii + 1)]:
            for value in line.split():
                numk += 1
                ffield[4][((ii + 1))][numk] = line.split()[line.split().index(value)]
                if numk > 2:
                    id_flag = "4-" + str(ii + 1) + '-' + str(numk - 2)
                    ff_value = float(line.split()[line.split().index(value)])
                    ffi[id_flag] = ff_value

    # sec 5
    numi = 0
    numj = 1
    numk = 0
    for ii in range(numsec5):
        ffield[5][(ii + 1)] = {}
        numk = 0
        for line in fflines[(sec5ids + ii): (sec5ids + ii + 1)]:
            for value in line.split():
                numk += 1
                ffield[5][((ii + 1))][numk] = line.split()[line.split().index(value)]
                if numk > 3:
                    id_flag = "5-" + str(ii + 1) + '-' + str(numk - 3)
                    ff_value = float(line.split()[line.split().index(value)])
                    ffi[id_flag] = ff_value

    # sec 6
    numi = 0
    numj = 1
    numk = 0
    for ii in range(numsec6):
        ffield[6][(ii + 1)] = {}
        numk = 0
        for line in fflines[(sec6ids + ii): (sec6ids + ii + 1)]:
            for value in line.split():
                numk += 1
                ffield[6][((ii + 1))][numk] = line.split()[line.split().index(value)]
                if numk > 4:
                    id_flag = "6-" + str(ii + 1) + '-' + str(numk - 4)
                    ff_value = float(line.split()[line.split().index(value)])
                    ffi[id_flag] = ff_value

    # sec 7
    numi = 0
    numj = 1
    numk = 0
    if numsec7 != 0:
        for ii in range(numsec7):
            ffield[7][(ii + 1)] = {}
            numk = 0
            for line in fflines[(sec7ids + ii): (sec7ids + ii + 1)]:
                for value in line.split():
                    numk += 1
                    ffield[7][((ii + 1))][numk] = line.split()[line.split().index(value)]
                    if numk > 3:
                        id_flag = "7-" + str(ii + 1) + '-' + str(numk - 3)
                        ff_value = float(line.split()[line.split().index(value)])
                        ffi[id_flag] = ff_value
    #
    ff_value_true = {}
    atom_type_id_dir = {}
    for id2 , id2_v in ffield[2].items():
        atom_id_num = id2
        atom_id_type = id2_v[1]
        atom_type_id_dir[str(atom_id_num)] = atom_id_type
    print(atom_type_id_dir)
    atom_type_id_dir["0"] = '*'
    for id1 ,sec_v in ffield.items():
        for id2 ,hang_v in sec_v.items():
            if id1 == 1:
                id2_flag = "0"
                id2_flag_id = ''
            elif id1 == 2:
                id2_flag = ffield[id1][id2][1]
                id2_flag_id = ''
            elif id1 == 3 or id1 == 4:
                id2_flag = atom_type_id_dir[ffield[id1][id2][1]] + "-" + atom_type_id_dir[ffield[id1][id2][2]]
                id2_flag_id = ffield[id1][id2][1] + "-" + ffield[id1][id2][2]
            elif id1 == 5 or id1 == 7:
                id2_flag = atom_type_id_dir[ffield[id1][id2][1]] + "-" + atom_type_id_dir[ffield[id1][id2][2]] \
                           +"-" + atom_type_id_dir[ffield[id1][id2][3]]
                id2_flag_id = ffield[id1][id2][1] + "-" + ffield[id1][id2][2] + "-" + \
                           ffield[id1][id2][3]
            elif id1 == 6 :
                id2_flag = atom_type_id_dir[ffield[id1][id2][1]] + "-" + atom_type_id_dir[ffield[id1][id2][2]]\
                           +"-" + atom_type_id_dir[ffield[id1][id2][3]] +"-" + atom_type_id_dir[ffield[id1][id2][4]]
                id2_flag_id = ffield[id1][id2][1] + "-" + ffield[id1][id2][2]\
                           +"-" + ffield[id1][id2][3] +"-" + ffield[id1][id2][4]
            for id3 , value in hang_v.items():
                if id1 == 1:
                    id_flag =  str(id1) + "-0-" + str(id3)
                elif id1 == 2:
                    id_flag =  str(id1) + "-" + str(id2) + '-' + str(id3 )
                elif id1 == 3 or id1 == 4:
                    id_flag = str(id1) + "-" + str(id2 ) + '-' + str(id3 )
                elif id1 == 5 or id1 == 7:
                    id_flag =  str(id1) + "-" + str(id2 ) + '-' + str(id3)
                elif id1 == 6:
                    id_flag =  str(id1) + "-" + str(id2 ) + '-' + str(id3 )
                id_flag_true = str(id1) + "_" + id2_flag + '  ' + id2_flag_id
                ff_value_true[id_flag] = id_flag_true

    params_file = path + os.sep + "params"
    with open(params_file, 'w') as w:
        w.write("")
    id_flag_1 = []
    id_flag_1.append("0")
    id_flag_2 = []
    id_flag_2.append("0")

    for i in pamlines:
        perpam = i.split()
        if len(perpam) > 0 and i.find("#") != 0:
            id1 = perpam[0]
            id2 = perpam[1]
            id3 = perpam[2]
            id_flag = id1 + '-' + id2 + '-' + id3
            v_a = float(params_start[id_flag]["a"])
            v_b = float(params_start[id_flag]["b"])
            v_v = float(ffi[id_flag])
            v_scale = v_b -v_a
            if v_scale == 0:
                print("params.scale  " + id_flag)
                sys.exit()
            # 1
            params_a = (v_v - scale * v_scale)
            params_b = (v_v + scale * v_scale)
            # 2
            if abs(v_v) != 0 and scale != 1:
                derta = abs(v_v) * scale
                params_a = (v_v - derta)
                params_b = (v_v + derta)
            #
            params_a = float("%.4f" % (max(params_a, v_a)))
            params_b = float("%.4f" % (min(params_b, v_b)))

            id_flag_1.append(id1)
            id_flag_2.append(id2)
            step = "%.3f" % ((params_b - params_a) / 20)
            if id_flag_1[-2] != id1:
                with open(params_file, 'a') as w:
                    w.write("#########\n" + "# "+ ff_value_true[id_flag]  + "\n" \
                            + str(id1).ljust(3) + str(id2).ljust(3) + str(id3).ljust(5) + str(v_v).ljust(10) + \
                            str(step).center(10) + str(params_a).ljust(10) + str(params_b).ljust(10) + '\n')
            else:
                if id_flag_2[-2] != id2:
                    with open(params_file, 'a') as w:
                        w.write("#\n" +"# "+ ff_value_true[id_flag]  + "\n" \
                                + str(id1).ljust(3) + str(id2).ljust(3) + str(id3).ljust(5) + str(v_v).ljust(10) + \
                                str(step).center(10) + str(params_a).ljust(10) + str(params_b).ljust(10) + '\n')
                else:
                    with open(params_file, 'a') as w:
                        w.write(str(id1).ljust(3) + str(id2).ljust(3) + str(id3).ljust(5) + str(v_v).ljust(10) + \
                                str(step).center(10) + str(params_a).ljust(10) + str(params_b).ljust(10) + '\n')
    return ff_value_true


