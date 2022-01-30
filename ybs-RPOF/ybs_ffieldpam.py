def ffpam(name,x,outfilename='outdata'):
    import os, sys
    import string
    from typing import List

    fname = name
    a = x
    outdir = outfilename
    dir_path = os.getcwd()
    if not os.path.exists(os.curdir + os.sep + outdir + os.sep + fname):
        ff = open(dir_path + os.sep + fname,'r')
    else:
        ff = open(dir_path + os.sep + outdir + os.sep + fname, 'r')
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
    for line in fflines:
        spline = line.split()
        if ((('!' in line) and type(eval(line.split()[0])) == float)):
            numj += 1
            ffield[1][0][numj] = spline[0]

    #sec2

    numi = 0
    numj = 1
    numk = 0
    for ii in range(numsec2):
        ffield[2][(ii)+1] = {}
        numk = 0
        for line in fflines[(sec2ids + ii * (4+lg) ): ((sec2ids + ii * (4+lg)) + (4+lg))]:
            for value in line.split():
                numk += 1
                ffield[2][(ii+1)][numk] = line.split()[line.split().index(value)]

    #sec3

    numi = 0
    numj = 1
    numk = 0
    for ii in range(numsec3):
        ffield[3][(ii+1)] = {}
        numk = 0
        for line in fflines[(sec3ids+ii*2) : ((sec3ids+ii*2)+2)]:
            for value in line.split():
                numk += 1
                ffield[3][((ii+1))][numk] = line.split()[line.split().index(value)]

    # sec 4
    numi = 0
    numj = 1
    numk = 0
    for ii in range(numsec4):
        ffield[4][(ii+1)] = {}
        numk = 0
        for line in fflines[(sec4ids+ii) : (sec4ids+ii +1)]:
            for value in line.split():
                numk += 1
                ffield[4][((ii+1))][numk] = line.split()[line.split().index(value)]

    # sec 5
    numi = 0
    numj = 1
    numk = 0
    for ii in range(numsec5):
        ffield[5][(ii+1)] = {}
        numk = 0
        for line in fflines[(sec5ids+ii) : (sec5ids+ii +1)]:
            for value in line.split():
                numk += 1
                ffield[5][((ii+1))][numk] = line.split()[line.split().index(value)]

    # sec 6
    numi = 0
    numj = 1
    numk = 0
    for ii in range(numsec6):
        ffield[6][(ii+1)] = {}
        numk = 0
        for line in fflines[(sec6ids+ii) : (sec6ids+ii +1)]:
            for value in line.split():
                numk += 1
                ffield[6][((ii+1))][numk] = line.split()[line.split().index(value)]

    # sec 7
    numi = 0
    numj = 1
    numk = 0
    if numsec7 != 0:
        for ii in range(numsec7):
            ffield[7][(ii+1)] = {}
            numk = 0
            for line in fflines[(sec7ids+ii) : (sec7ids+ii +1)]:
                for value in line.split():
                    numk += 1
                    ffield[7][((ii+1))][numk] = line.split()[line.split().index(value)]


    if (a[0] == 2):
        a[2] += 1
    elif (a[0] == 3) or (a[0] == 4):
        a[2] += 2
    elif (a[0] == 5) or (a[0] == 7):
        a[2] += 3
    elif (a[0] == 6):
        a[2] += 4
    try:
        return_v = float(ffield[a[0]][a[1]][a[2]])
    except:
        return_v = 999999
    return return_v

