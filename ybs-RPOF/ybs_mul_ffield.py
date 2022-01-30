import os, sys, time
import string
from typing import List
def mul_ffield(x_pam,pso_ff_name,outfilename='outdata'):


    x_flag_value = x_pam
    ff_source_name = "ffield.reax.start"
    outdir = outfilename
    ff_pso_name = pso_ff_name
    dir_path = os.getcwd()
    try:
        ff = open(dir_path + os.sep + outdir + os.sep + ff_source_name, 'r')
    except:
        sys.exit()
    fflines = ff.readlines()
    ffield = {}

    #

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
        if ((('!' in line) and type(eval(line.split()[0])) == float)) or sec2tid > fflines.index(line) >= sec1ids:
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
    if type(x_flag_value) is list:
        for kk in x_flag_value:
            id1 = int(kk[0])
            id2 = int(kk[1])
            id3 = int(kk[2])
            pam_value = kk[3]
            pam_id1 = id1
            pam_id2 = id2
            if id1 == 1:
                pam_id3 = id3
            elif id1 == 2:
                pam_id3 = id3 + 1
            elif (id1 == 3) or (id1 == 4):
                pam_id3 = id3 + 2
            elif (id1 == 5) or (id1 == 7):
                pam_id3 = id3 + 3
            elif (id1 == 6):
                pam_id3 = id3 + 4
            ffield[pam_id1][pam_id2][pam_id3] = "%.4f" % float(pam_value)
    else:
        for kk,vv in x_flag_value.items():
            pam_flag = kk
            pam_value = vv
            id1 = int(pam_flag.split("-")[0])
            id2 = int(pam_flag.split("-")[1])
            id3 = int(pam_flag.split("-")[2])
            pam_id1 = id1
            pam_id2 = id2
            if id1 == 1:
                pam_id3 = id3
            elif id1 == 2:
                pam_id3 = id3 + 1
            elif (id1 == 3) or (id1 == 4):
                pam_id3 = id3 + 2
            elif (id1 == 5) or (id1 == 7):
                pam_id3 = id3 + 3
            elif (id1 == 6):
                pam_id3 = id3 + 4
            ffield[pam_id1][pam_id2][pam_id3] = "%.4f" % float(pam_value)

    ffnew = open(dir_path +os.sep+ outdir +os.sep + ff_pso_name,'w')
    ffnew.write("NEWffield " + time.strftime("%Y-%m-%d %I:%M:%S") + '\n')
    ffnew.close()
    ffnew = open(dir_path +os.sep+ outdir +os.sep + ff_pso_name,'a')

    sec1_biaoshi = {1: "p(boc1)", 2: "p(boc2)", 3: "p(coa2)", 4: "p(trip4)", 5: "p(trip3)", 6: "kc2", 7: "p(ovun6)",
                    8: "p(trip2)", 9: "p(ovun7)", 10: "p(ovun8)", 11: "p(trip1)", 12: "Lower Taper-radius (swa)",
                    13: "Upper Taper-radius (swb)", 14: "not used", 15: "p(val7)", 16: "p(lp1)", 17: "p(val9)",
                    18: "p(val10)", 19: "not used", 20: "p(pen2)", 21: "p(pen3)", 22: "p(pen4)", 23: "not used",
                    24: "p(tor2)", 25: "p(tor3)", 26: "p(tor4)", 27: "not used", 28: "p(cot2)", 29: "p(vdW1)",
                    30: "Cutoff for bond order*100 (cutoff)", 31: "p(coa4)", 32: "p(ovun4)", 33: "p(ovun3)",
                    34: "p(val8)", 35: "not used", 36: "not used", 37: "not used", 38: "not used", 39: "p(coa3)"}

    #sec 1
    i = 0
    j = 0
    k = 0
    ffnew.write(sectitle[0] )
    for num in range(numsec1):
        i = 1
        k += 1
        ffnew.write(ffield[i][0][k].rjust(12) + "  ! " + sec1_biaoshi[k] + '\n')

    #sec 2
    i = 0
    j = 0
    k = 0
    for line in fflines[sec2tid : atomid[0]]:
        ffnew.write(line)
    for num in range(numsec2):
        i = 2
        j += 1
        ffnew.write(ffield[i][j][1].rjust(4) +ffield[i][j][2].rjust(10) + ffield[i][j][3].rjust(10)+ffield[i][j][4].rjust(10)+ffield[i][j][5].rjust(10)+ffield[i][j][6].rjust(10)+ffield[i][j][7].rjust(10)+ffield[i][j][8].rjust(10)+ffield[i][j][9].rjust(10)+'\n')
        ffnew.write(ffield[i][j][10].rjust(14) +ffield[i][j][11].rjust(10) + ffield[i][j][12].rjust(10)+ffield[i][j][13].rjust(10)+ffield[i][j][14].rjust(10)+ffield[i][j][15].rjust(10)+ffield[i][j][16].rjust(10)+ffield[i][j][17].rjust(10)+'\n')
        ffnew.write(ffield[i][j][18].rjust(14) +ffield[i][j][19].rjust(10) + ffield[i][j][20].rjust(10)+ffield[i][j][21].rjust(10)+ffield[i][j][22].rjust(10)+ffield[i][j][23].rjust(10)+ffield[i][j][24].rjust(10)+ffield[i][j][25].rjust(10)+'\n')
        ffnew.write(ffield[i][j][26].rjust(14) +ffield[i][j][27].rjust(10) + ffield[i][j][28].rjust(10)+ffield[i][j][29].rjust(10)+ffield[i][j][30].rjust(10)+ffield[i][j][31].rjust(10)+ffield[i][j][32].rjust(10)+ffield[i][j][33].rjust(10)+'\n')
        if lg == 1:
            ffnew.write(ffield[i][j][34].rjust(14) +ffield[i][j][35].rjust(10) +'\n')

    # sec 3
    i = 0
    j = 0
    k = 0
    for line in fflines[sec3tid : sec3ids]:
        ffnew.write(line)
    for num in range(numsec3):
        i = 3
        j +=1
        ffnew.write(ffield[i][j][1].rjust(2) +ffield[i][j][2].rjust(3) + ffield[i][j][3].rjust(9)+ffield[i][j][4].rjust(10)+ffield[i][j][5].rjust(10)+ffield[i][j][6].rjust(10)+ffield[i][j][7].rjust(10)+ffield[i][j][8].rjust(10)+ffield[i][j][9].rjust(10)+ffield[i][j][10].rjust(10)+'\n')
        ffnew.write(ffield[i][j][11].rjust(14) +ffield[i][j][12].rjust(10) + ffield[i][j][13].rjust(10)+ffield[i][j][14].rjust(10)+ffield[i][j][15].rjust(10)+ffield[i][j][16].rjust(10)+ffield[i][j][17].rjust(10)+ffield[i][j][18].rjust(10)+'\n')

    # sec 4
    i = 0
    j = 0
    k = 0
    ffnew.write(fflines[sec4tid])
    for num in range(numsec4):
        i = 4
        j +=1
        if lg == 1:
            ffnew.write(ffield[i][j][1].rjust(2) +ffield[i][j][2].rjust(3) + ffield[i][j][3].rjust(9)+ffield[i][j][4].rjust(10)+ffield[i][j][5].rjust(10)+ffield[i][j][6].rjust(10)+ffield[i][j][7].rjust(10)+ffield[i][j][8].rjust(10)+ffield[i][j][9].rjust(11)+'\n')
        else:
            ffnew.write(ffield[i][j][1].rjust(2) +ffield[i][j][2].rjust(3) + ffield[i][j][3].rjust(9)+ffield[i][j][4].rjust(10)+ffield[i][j][5].rjust(10)+ffield[i][j][6].rjust(10)+ffield[i][j][7].rjust(10)+ffield[i][j][8].rjust(10)+'\n')

    # sec 5
    i = 0
    j = 0
    k = 0
    ffnew.write(fflines[sec5tid])
    for num in range(numsec5):
        i = 5
        j +=1
        ffnew.write(ffield[i][j][1].rjust(2) +ffield[i][j][2].rjust(3) + ffield[i][j][3].rjust(3)+ffield[i][j][4].rjust(10)+ffield[i][j][5].rjust(10)+ffield[i][j][6].rjust(10)+ffield[i][j][7].rjust(10)+ffield[i][j][8].rjust(10)+ffield[i][j][9].rjust(10)+ffield[i][j][10].rjust(10)+'\n')

    # sec 6
    i = 0
    j = 0
    k = 0
    ffnew.write(fflines[sec6tid])
    for num in range(numsec6):
        i = 6
        j +=1
        ffnew.write(ffield[i][j][1].rjust(2) +ffield[i][j][2].rjust(3) + ffield[i][j][3].rjust(3)+ffield[i][j][4].rjust(3)+ffield[i][j][5].rjust(10)+ffield[i][j][6].rjust(10)+ffield[i][j][7].rjust(10)+ffield[i][j][8].rjust(10)+ffield[i][j][9].rjust(10)+ffield[i][j][10].rjust(10)+ffield[i][j][10].rjust(10)+'\n')

    # sec 7
    i = 0
    j = 0
    k = 0
    ffnew.write(fflines[sec7tid])
    if numsec7 != 0:
        for num in range(numsec7):
            i = 7
            j +=1
            ffnew.write(ffield[i][j][1].rjust(2) +ffield[i][j][2].rjust(3) + ffield[i][j][3].rjust(3)+ffield[i][j][4].rjust(10)+ffield[i][j][5].rjust(10)+ffield[i][j][6].rjust(10)+ffield[i][j][7].rjust(10)+'\n')
    ffnew.close()
    ff.close()
    return ff_pso_name

