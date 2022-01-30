import string
from ybs_trainset import trainset as ts
def optgeo():
    tset, charge, geometry, energy, cellpam, heat = ts('trainset')
    optname = []
    optcul = {}
    energygeoc = []

    chargename =[]

    for num ,v in tset['CHARGE'].items():
        chargename.append(tset['CHARGE'][num]['geo'])
    [optname.append(i) for i in chargename if not i in optname ]
    geometryname = []

    for num,v in tset['GEOMETRY'].items():
        geometryname.append(tset['GEOMETRY'][num]['geo'])
    [optname.append(i) for i in geometryname if not i in optname ]

    energyname = []
    for num,v in tset['ENERGY'].items():
        if tset['ENERGY'][num]['numgeo'] == 2:
            energyname.append(tset['ENERGY'][num]['geo2'])
        else:
            continue
    [optname.append(i) for i in energyname if not i in optname ]

    q = []
    x = []
    reax = 0.0
    for i in optname:
        optcul[i] = {}
        optcul[i]['q']= q
        optcul[i]['x']= x
        optcul[i]['reax']= reax

    energygeo = []
    for key,value in tset['ENERGY'].items():
        if value['numgeo'] == 1:
            energygeo.append(value['geo1'])
        if value['numgeo'] == 2:
            energygeo.append(value['geo1'])
            energygeo.append(value['geo2'])
        if value['numgeo'] == 3:
            energygeo.append(value['geo1'])
            energygeo.append(value['geo2'])
            energygeo.append(value['geo3'])
        if value['numgeo'] == 4:
            energygeo.append(value['geo1'])
            energygeo.append(value['geo2'])
            energygeo.append(value['geo3'])
            energygeo.append(value['geo4'])
        if value['numgeo'] == 5:
            energygeo.append(value['geo1'])
            energygeo.append(value['geo2'])
            energygeo.append(value['geo3'])
            energygeo.append(value['geo4'])
            energygeo.append(value['geo5'])
    energygeoa=[]
    [energygeoa.append(i) for i in energygeo if not i in optname ]
    [energygeoc.append(i) for i in energygeoa if not i in energygeoc ]
    #v=[[0 for i in range(3)] for i in range(3)]
    #optcul['H3SiCl']['x']=v
    return(optcul,energygeoc)

