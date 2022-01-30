import os, sys, re
import string
from ybs_bgf import bgfdata
def trainset(setname):
    getname = setname
    dir_path = os.getcwd()
    sf = open(dir_path +os.sep +  getname,'r')
    setlines = sf.readlines()
    numcharge = 0
    numgeometry = 0
    numenergy = 0
    numcellpam = 0
    numheat = 0
    numreaction_energy = 0
    trainset = {}
    charge = {}
    geometry = {}
    energy ={}
    reaction_energy = {}
    cellpam = {}
    heat = {}

    chargesid = 0
    chargeeid = 0
    geometrysid = 0
    geometryeid = 0
    energysid = 0
    energyeid = 0
    cellpamsid = 0
    cellpameid = 0
    heatsid = 0
    heateid = 0
    reaction_energysid = 0
    reaction_energyeid = 0
    for line in setlines:
        if line.find('CHARGE') == 0:
            name=line.split()[0]
            trainset[name] = {}
            chargesid = setlines.index(line) + 1
        if line.find('ENDCHARGE') == 0:
            chargeeid = setlines.index(line)

        if line.find('GEOMETRY') == 0:
            name=line.split()[0]
            trainset[name] = {}
            geometrysid = setlines.index(line) + 1
        if line.find('ENDGEOMETRY') == 0:
            geometryeid = setlines.index(line)

        if line.find('ENERGY') == 0:
            name=line.split()[0]
            trainset[name] = {}
            energysid = setlines.index(line) + 1
        if line.find('ENDENERGY') == 0:
            energyeid = setlines.index(line)
        # reaction_energy
        if line.find('REACTION_ENERGY') == 0:
            name=line.split()[0]
            trainset[name] = {}
            reaction_energysid = setlines.index(line) + 1
        if line.find('ENDREACTION_ENERGY') == 0:
            reaction_energyeid = setlines.index(line)

        if line.find('CELL_PARAMETERS') == 0:
            name=line.split()[0]
            trainset[name] = {}
            cellpamsid = setlines.index(line) + 1
        if line.find('ENDCELL_PARAMETERS') == 0:
            cellpameid = setlines.index(line)

        if line.find('HEATFO') == 0:
            name=line.split()[0]
            trainset[name] = {}
            heatsid = setlines.index(line) + 1
        if line.find('ENDHEATFO') == 0:
            heateid = setlines.index(line)

    # 'CHARGE
    charge['geo'] = []
    for line in setlines[chargesid:chargeeid]:
        spline = line.split()
        if line.find("#") != 0 and len(spline) > 0 :
            numcharge += 1
            trainset['CHARGE'][numcharge]={}
            trainset['CHARGE'][numcharge]['geo']= spline[0]
            trainset['CHARGE'][numcharge]['weight']= float(spline[1])
            trainset['CHARGE'][numcharge]['atom']= int(spline[2])
            trainset['CHARGE'][numcharge]['value']= float(spline[-1])
            charge['geo'].append(spline[0])
    charge['num'] = numcharge

    #CELL_PARAMETERS
    cellpam['geo'] = []
    for line in setlines[cellpamsid:cellpameid]:
        spline = line.split()
        if line.find("#") != 0 and len(spline) > 0:
            numcellpam += 1
            trainset['CELL_PARAMETERS'][numcellpam] = {}
            trainset['CELL_PARAMETERS'][numcellpam]['geo'] = spline[0]
            trainset['CELL_PARAMETERS'][numcellpam]['weight'] = float(spline[1])
            trainset['CELL_PARAMETERS'][numcellpam]['type'] = spline[2]
            trainset['CELL_PARAMETERS'][numcellpam]['value'] = float(spline[-1])
            cellpam['geo'].append(spline[0])
    cellpam['num'] = numcellpam

    #HEATFO
    heat['geo'] = []
    for line in setlines[heatsid:heateid]:
        spline = line.split()
        if line.find("#") != 0 and len(spline) > 0:
            numheat += 1
            trainset['HEATFO'][numheat] = {}
            trainset['HEATFO'][numheat]['geo'] = spline[0]
            trainset['HEATFO'][numheat]['weight'] = float(spline[1])
            trainset['HEATFO'][numheat]['value'] = float(spline[-1])
            heat['geo'].append(spline[0])
    heat['num'] = numheat

    # 'GEOMETRY'
    geometry['geo'] = []
    for line in setlines[geometrysid:geometryeid]:
        spline = line.split()
        if line.find("#") != 0 and len(spline) > 0:
            numgeometry += 1
            trainset['GEOMETRY'][numgeometry]={}
            trainset['GEOMETRY'][numgeometry]['geo']= spline[0]
            trainset['GEOMETRY'][numgeometry]['weight']= float(spline[1])
            trainset['GEOMETRY'][numgeometry]['value']= float(spline[-1])
            trainset['GEOMETRY'][numgeometry]['atom1']= int(spline[2])
            trainset['GEOMETRY'][numgeometry]['atom2']= int(spline[3])
            geometry['geo'].append(spline[0])
            if len(spline) == 5:
                trainset['GEOMETRY'][numgeometry]['numatom']= 2
                trainset['GEOMETRY'][numgeometry]['type']= 'bond'

            if len(spline) == 6:
                trainset['GEOMETRY'][numgeometry]['numatom']= 3
                trainset['GEOMETRY'][numgeometry]['type']= 'angle'
                trainset['GEOMETRY'][numgeometry]['atom3']= int(spline[4])
            if len(spline) == 7:
                trainset['GEOMETRY'][numgeometry]['numatom']= 4
                trainset['GEOMETRY'][numgeometry]['type']= 'torsion'
                trainset['GEOMETRY'][numgeometry]['atom3']= int(spline[4])
                trainset['GEOMETRY'][numgeometry]['atom4']= int(spline[5])
    geometry['num'] = numgeometry

    # 'ENERGY'
    energy['geo'] = []
    for line in setlines[energysid:energyeid]:
        spline_ = line.split()
        if line.find("#") != 0 and len(spline_) > 0:
            line = line[:-2]
            numenergy += 1
            spline = re.split('[ |/]+',line)
            numsp = int((len(spline) - 2)/3)
            trainset['ENERGY'][numenergy]={}
            trainset['ENERGY'][numenergy]['weight']= float(spline[0])
            trainset['ENERGY'][numenergy]['value']= float(spline[-1])
            trainset['ENERGY'][numenergy]['geo1']= spline[2]
            trainset['ENERGY'][numenergy]['geo1+']= spline[1]
            trainset['ENERGY'][numenergy]['geo1/']= float(spline[3])
            energy['geo'].append(spline[2])
            if numsp == 1:
                trainset['ENERGY'][numenergy]['numgeo']= 1
            if numsp == 2:
                trainset['ENERGY'][numenergy]['numgeo']= 2
                trainset['ENERGY'][numenergy]['geo2+']= spline[4]
                trainset['ENERGY'][numenergy]['geo2']= spline[5]
                trainset['ENERGY'][numenergy]['geo2/']= float(spline[6])
                energy['geo'].append(spline[5])
            if numsp == 3:
                trainset['ENERGY'][numenergy]['numgeo']= 3
                trainset['ENERGY'][numenergy]['geo2+']= spline[4]
                trainset['ENERGY'][numenergy]['geo2']= spline[5]
                trainset['ENERGY'][numenergy]['geo2/']= float(spline[6])
                trainset['ENERGY'][numenergy]['geo3+']= spline[7]
                trainset['ENERGY'][numenergy]['geo3']= spline[8]
                trainset['ENERGY'][numenergy]['geo3/']= float(spline[9])
                energy['geo'].append(spline[5])
                energy['geo'].append(spline[8])
            if numsp == 4:
                trainset['ENERGY'][numenergy]['numgeo']= 4
                trainset['ENERGY'][numenergy]['geo2+']= spline[4]
                trainset['ENERGY'][numenergy]['geo2']= spline[5]
                trainset['ENERGY'][numenergy]['geo2/']= float(spline[6])
                trainset['ENERGY'][numenergy]['geo3+']= spline[7]
                trainset['ENERGY'][numenergy]['geo3']= spline[8]
                trainset['ENERGY'][numenergy]['geo3/']= float(spline[9])
                trainset['ENERGY'][numenergy]['geo4+']= spline[10]
                trainset['ENERGY'][numenergy]['geo4']= spline[11]
                trainset['ENERGY'][numenergy]['geo4/']= float(spline[12])
                energy['geo'].append(spline[5])
                energy['geo'].append(spline[8])
                energy['geo'].append(spline[11])
            if numsp == 5:
                trainset['ENERGY'][numenergy]['numgeo']= 5
                trainset['ENERGY'][numenergy]['geo2+']= spline[4]
                trainset['ENERGY'][numenergy]['geo2']= spline[5]
                trainset['ENERGY'][numenergy]['geo2/']= float(spline[6])
                trainset['ENERGY'][numenergy]['geo3+']= spline[7]
                trainset['ENERGY'][numenergy]['geo3']= spline[8]
                trainset['ENERGY'][numenergy]['geo3/']= float(spline[9])
                trainset['ENERGY'][numenergy]['geo4+']= spline[10]
                trainset['ENERGY'][numenergy]['geo4']= spline[11]
                trainset['ENERGY'][numenergy]['geo4/']= float(spline[12])
                trainset['ENERGY'][numenergy]['geo5+']= spline[13]
                trainset['ENERGY'][numenergy]['geo5']= spline[14]
                trainset['ENERGY'][numenergy]['geo5/']= float(spline[15])
                energy['geo'].append(spline[5])
                energy['geo'].append(spline[8])
                energy['geo'].append(spline[11])
                energy['geo'].append(spline[14])
    energy['num'] = numenergy

    # reaction_energy
    reaction_energy['geo'] = []
    for line in setlines[reaction_energysid:reaction_energyeid]:
        spline_ = line.split()
        if line.find("#") != 0 and len(spline_) > 0:
            numreaction_energy += 1
            spline = re.split('[ |/]+',line)
            numsp = int((len(spline) - 2)/3)

            trainset['REACTION_ENERGY'][numreaction_energy]={}
            trainset['REACTION_ENERGY'][numreaction_energy]['weight']= float(spline[0])
            trainset['REACTION_ENERGY'][numreaction_energy]['value']= float(spline[-1])
            trainset['REACTION_ENERGY'][numreaction_energy]['geo1']= spline[2]
            trainset['REACTION_ENERGY'][numreaction_energy]['geo1+']= spline[1]
            trainset['REACTION_ENERGY'][numreaction_energy]['geo1/']= float(spline[3])
            reaction_energy['geo'].append(spline[2])
            if numsp == 1:
                trainset['REACTION_ENERGY'][numreaction_energy]['numgeo']= 1
            if numsp == 2:
                trainset['REACTION_ENERGY'][numreaction_energy]['numgeo']= 2
                trainset['REACTION_ENERGY'][numreaction_energy]['geo2+']= spline[4]
                trainset['REACTION_ENERGY'][numreaction_energy]['geo2']= spline[5]
                trainset['REACTION_ENERGY'][numreaction_energy]['geo2/']= float(spline[6])
                reaction_energy['geo'].append(spline[5])
            if numsp == 3:
                trainset['REACTION_ENERGY'][numreaction_energy]['numgeo']= 3
                trainset['REACTION_ENERGY'][numreaction_energy]['geo2+']= spline[4]
                trainset['REACTION_ENERGY'][numreaction_energy]['geo2']= spline[5]
                trainset['REACTION_ENERGY'][numreaction_energy]['geo2/']= float(spline[6])
                trainset['REACTION_ENERGY'][numreaction_energy]['geo3+']= spline[7]
                trainset['REACTION_ENERGY'][numreaction_energy]['geo3']= spline[8]
                trainset['REACTION_ENERGY'][numreaction_energy]['geo3/']= float(spline[9])
                reaction_energy['geo'].append(spline[5])
                reaction_energy['geo'].append(spline[8])
            if numsp == 4:
                trainset['REACTION_ENERGY'][numreaction_energy]['numgeo']= 4
                trainset['REACTION_ENERGY'][numreaction_energy]['geo2+']= spline[4]
                trainset['REACTION_ENERGY'][numreaction_energy]['geo2']= spline[5]
                trainset['REACTION_ENERGY'][numreaction_energy]['geo2/']= float(spline[6])
                trainset['REACTION_ENERGY'][numreaction_energy]['geo3+']= spline[7]
                trainset['REACTION_ENERGY'][numreaction_energy]['geo3']= spline[8]
                trainset['REACTION_ENERGY'][numreaction_energy]['geo3/']= float(spline[9])
                trainset['REACTION_ENERGY'][numreaction_energy]['geo4+']= spline[10]
                trainset['REACTION_ENERGY'][numreaction_energy]['geo4']= spline[11]
                trainset['REACTION_ENERGY'][numreaction_energy]['geo4/']= float(spline[12])
                reaction_energy['geo'].append(spline[5])
                reaction_energy['geo'].append(spline[8])
                reaction_energy['geo'].append(spline[11])
            if numsp == 5:
                trainset['REACTION_ENERGY'][numreaction_energy]['numgeo']= 5
                trainset['REACTION_ENERGY'][numreaction_energy]['geo2+']= spline[4]
                trainset['REACTION_ENERGY'][numreaction_energy]['geo2']= spline[5]
                trainset['REACTION_ENERGY'][numreaction_energy]['geo2/']= float(spline[6])
                trainset['REACTION_ENERGY'][numreaction_energy]['geo3+']= spline[7]
                trainset['REACTION_ENERGY'][numreaction_energy]['geo3']= spline[8]
                trainset['REACTION_ENERGY'][numreaction_energy]['geo3/']= float(spline[9])
                trainset['REACTION_ENERGY'][numreaction_energy]['geo4+']= spline[10]
                trainset['REACTION_ENERGY'][numreaction_energy]['geo4']= spline[11]
                trainset['REACTION_ENERGY'][numreaction_energy]['geo4/']= float(spline[12])
                trainset['REACTION_ENERGY'][numreaction_energy]['geo5+']= spline[13]
                trainset['REACTION_ENERGY'][numreaction_energy]['geo5']= spline[14]
                trainset['REACTION_ENERGY'][numreaction_energy]['geo5/']= float(spline[15])
                reaction_energy['geo'].append(spline[5])
                reaction_energy['geo'].append(spline[8])
                reaction_energy['geo'].append(spline[11])
                reaction_energy['geo'].append(spline[14])
    reaction_energy['num'] = numreaction_energy
    ##
    sf.close()

    chageo = []
    gemgeo = []
    enegeo = []
    #
    reacenegeo = []
    cellpamgeo = []
    heatgeo = []
    nominigeo = []
    minigeo = []
    [chageo.append(i) for i in charge['geo'] if not i in chageo]
    [enegeo.append(i) for i in energy['geo'] if not i in enegeo]
    # reaction_energy
    [reacenegeo.append(i) for i in reaction_energy['geo'] if not i in reacenegeo]
    [gemgeo.append(i) for i in geometry['geo'] if not i in gemgeo]
    [cellpamgeo.append(i) for i in cellpam['geo'] if not i in cellpamgeo]
    [heatgeo.append(i) for i in heat['geo'] if not i in heatgeo]
    #
    sumgeo = chageo + enegeo #+ heatgeo
    sumgeo1 = gemgeo + cellpamgeo + reacenegeo + heatgeo
    # nomini
    [nominigeo.append(i) for i in sumgeo if not i in nominigeo]
    # mini
    [minigeo.append(i) for i in sumgeo1 if not i in minigeo]
    #
    trainset_geo_list = []
    zonghe = sumgeo + sumgeo1
    [trainset_geo_list.append(i) for i in zonghe if not i in trainset_geo_list]
    bgf_geo_list = bgfdata("geofile")[0]
    for geo_name in trainset_geo_list:
        try:
            bgf_geo_list[geo_name]
        except:
            print("\n\nerror：  " + geo_name + "     \n\n")
    return(trainset, charge, geometry, energy,reaction_energy, cellpam, heat, nominigeo, minigeo)