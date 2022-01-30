import string
import threading
import multiprocessing
import math
from lammps import lammps
from ybs_ffieldheat import heat
def fderta(ffi,ts,gnt,fsum):
    #
    ffield = ffi
    tset = ts
    geontype = gnt
    fnomini = fsum[0]
    fmini   = fsum[1]
    #charge
    dertacharge = {}
    formcharge = {}
    if 'CHARGE' in tset:
        for key,value in tset['CHARGE'].items():
            num = key
            dertacharge[num] ={}
            dertacharge[num]['weight'] = value['weight']
            dertacharge[num]['derta'] = math.fabs(value['value']- (fnomini[value['geo']]['q'][value['atom']-1]))
            dertacharge[num]['QM-value'] = value['value']
            formcharge[num] = str(value['weight']).ljust(10) + value['geo'].ljust(30) + \
                              str(value['atom']).ljust(10) + str('%.4f'% value['value']).rjust(10) +\
                              str('%.4f'%(fnomini[value['geo']]['q'][value['atom']-1])).rjust(10) +\
                              str('%.4f'% dertacharge[num]['derta']).rjust(10)

    #CELL_PARAMETERS
    dertacellpam = {}
    formcellpam = {}
    if 'CELL_PARAMETERS' in tset:
        for key,value in tset['CELL_PARAMETERS'].items():
            num = key
            dertacellpam[num] ={}
            dertacellpam[num]['weight'] = value['weight']
            dertacellpam[num]['derta'] = math.fabs(value['value']- (fmini[value['geo']][value['type']]))
            dertacellpam[num]['QM-value'] = value['value']
            formcellpam[num] = str(value['weight']).ljust(10) + value['geo'].ljust(30) + \
                              str(value['type']).ljust(10) + str('%.4f'% value['value']).rjust(10) +\
                              str('%.4f'%(fmini[value['geo']][value['type']])).rjust(10) +\
                              str('%.4f'% dertacellpam[num]['derta']).rjust(10)

    # HEATFO
    dertaheat = {}
    formheat = {}
    heat_type = {"C": 183.8108, "H": 58.4369, "O": 69.2921, "N": 127.9672}
    #heat_type = {"C": 173.00, "H": 52.00, "O": 60.00, "N": 120.00}  #120
    if 'HEATFO' in tset:
        for key, value in tset['HEATFO'].items():
            num = key
            geo = value['geo']
            sumheat =0.0
            sumheat_gas = 0.0
            # geo - N2 -O2 -H2 - C
            enrgy_H2 = fmini["H2"]['energy'] / geontype["H2"]["H"]
            enrgy_O2 = fmini["O2"]['energy'] / geontype["O2"]["O"]
            enrgy_N2 = fmini["N2"]['energy'] / geontype["N2"]["N"]
            #enrgy_C = fmini["graphite-221"]['energy'] / geontype["graphite-221"]["C"]
            enrgy_C = -199.67333192597138
            heat_gas_type = {"C": enrgy_C , "H": enrgy_H2,"O":enrgy_O2 , "N":enrgy_N2}

            for kk, vv in geontype[geo].items():
                heat1 = heat_type[kk] #heat(ffield,kk)
                heat2 = heat1 * vv
                sumheat += heat2
                # geo - N2 -O2 -H2 - C
                heat_1 = heat_gas_type[kk]
                heat_2 = heat_1 * vv
                sumheat_gas += heat_2

            Hf = fmini[geo]['energy'] - sumheat_gas
            #Hf = fmini[geo]['energy'] + (4 * 8.3145 * 298.15 / 4185.8518208) + sumheat
            dertaheat[num] ={}
            dertaheat[num]['weight'] = value['weight']
            dertaheat[num]['derta'] = math.fabs(value['value'] - Hf)
            dertaheat[num]['QM-value'] = value['value']
            formheat[num] =   str(value['weight']).ljust(10) + value['geo'].ljust(30) +\
                              str('%.4f'% value['value']).rjust(10) +\
                              str('%.4f'%(Hf)).rjust(10) +\
                              str('%.4f'% dertaheat[num]['derta']).rjust(10)

    #geometry
    dertageometry = {}
    formgeometry = {}
    if 'GEOMETRY' in tset:
        from ybs_culgeo import bond as cb
        from ybs_culgeo import angle as ca
        from ybs_culgeo import torsion as ct

        for key,value in tset['GEOMETRY'].items():
            num = key
            atom1 = []
            atom2 = []
            atom3 = []
            atom4 = []
            bond = 0
            angle = 0
            torsion =0
            if value['type'] == 'bond':
                atom1 = fmini[value['geo']]['x'][value['atom1']-1]
                atom2 = fmini[value['geo']]['x'][value['atom2']-1]
                bond = cb(atom1,atom2)
                dertageometry[num] ={}
                dertageometry[num]['weight'] = value['weight']
                dertageometry[num]['derta'] = math.fabs(bond - value['value'])
                dertageometry[num]['QM-value'] = value['value']
                formgeometry[num] = str(value['weight']).ljust(10) + value['geo'].ljust(30) + \
                                    str(value['atom1']).ljust(3) +  str(value['atom2']).ljust(12) +\
                                    str('%.4f'% value['value']).rjust(10) +\
                                    str('%.4f'% bond).rjust(10) +\
                                    str('%.4f'% dertageometry[num]['derta']).rjust(10)

            elif value['type'] == 'angle':
                atom1 = fmini[value['geo']]['x'][value['atom1']-1]
                atom2 = fmini[value['geo']]['x'][value['atom2']-1]
                atom3 = fmini[value['geo']]['x'][value['atom3']-1]
                angle = ca(atom1,atom2,atom3)
                dertageometry[num] ={}
                dertageometry[num]['weight'] = value['weight']
                dertageometry[num]['derta'] = math.fabs(angle - value['value'])
                dertageometry[num]['QM-value'] = value['value']
                formgeometry[num] = str(value['weight']).ljust(10) + value['geo'].ljust(30) +  str(value['atom1']).ljust(3) +\
                                    str(value['atom2']).ljust(3) + str(value['atom3']).ljust(9) +\
                                    str('%.4f'% value['value']).rjust(10) +\
                                    str('%.4f'% angle).rjust(10) +\
                                    str('%.4f'% dertageometry[num]['derta']).rjust(10)

            elif value['type'] == 'torsion':
                atom1 = fmini[value['geo']]['x'][value['atom1']-1]
                atom2 = fmini[value['geo']]['x'][value['atom2']-1]
                atom3 = fmini[value['geo']]['x'][value['atom3']-1]
                atom4 = fmini[value['geo']]['x'][value['atom4']-1]
                torsion = ct(atom1,atom2,atom3,atom4)
                dertageometry[num] ={}
                dertageometry[num]['weight'] = value['weight']
                if (torsion * value['value'] < 0) and (math.fabs(torsion) > 100) and (math.fabs(value['value']) > 100):
                    dertageometry[num]['derta'] = math.fabs(math.fabs(torsion) - 180) + math.fabs(math.fabs(value['value']) - 180)
                elif (160 < math.fabs(value['value']) or math.fabs(value['value']) < 20) and math.fabs(torsion - value['value']) > 160:
                    dertageometry[num]['derta'] = 2 * math.fabs(math.fabs(torsion - value['value']) -180)
                else:
                    dertageometry[num]['derta']= math.fabs(torsion - value['value'])
                dertageometry[num]['QM-value'] = value['value']
                formgeometry[num] = str(value['weight']).ljust(10) + value['geo'].ljust(30) +  str(value['atom1']).ljust(3) +\
                                    str(value['atom2']).ljust(3) + str(value['atom3']).ljust(3) + str(value['atom4']).ljust(6) +\
                                    str('%.4f'% value['value']).rjust(10) +\
                                    str('%.4f'% torsion).rjust(10) +\
                                    str('%.4f'% dertageometry[num]['derta']).rjust(10)

    # energy
    dertaenergy = {}
    formenergy = {}
    if 'ENERGY' in tset:

        for key,value in tset['ENERGY'].items():
            num = key
            numg = value['numgeo']
            formula = []
            valenergy  = 0
            valenergy1 = 0
            valenergy2 = 0
            valenergy3 = 0
            valenergy4 = 0
            valenergy5 = 0
            if numg == 1:
                formula = value['geo1+'].ljust(2) + value['geo1'] + '/' + str(value['geo1/']).ljust(5)
                valenergy1 = fnomini[value['geo1']]['energy']
                valenergy = valenergy1
                dertaenergy[num] ={}
                dertaenergy[num]['weight'] = value['weight']
                dertaenergy[num]['derta'] = math.fabs(valenergy - value['value'])
                dertaenergy[num]['QM-value'] = value['value']
                formenergy[num] = str(value['weight']).ljust(10) + formula.ljust(60) +\
                                  str('%.4f'% value['value']).rjust(10) +\
                                  str('%.4f'% valenergy).rjust(10) +\
                                  str('%.4f'% dertaenergy[num]['derta']).rjust(10)

            if numg == 2:
                formula = (value['geo1+'].ljust(2) + value['geo1'] + '/' + str(value['geo1/']).ljust(5)) +\
                          (value['geo2+'].ljust(2) + value['geo2'] + '/' + str(value['geo2/']).ljust(5))

                valenergy1 = fnomini[value['geo1']]['energy']
                valenergy2 = fnomini[value['geo2']]['energy']
                valenergy = eval('%s%f/%f %s%f/%f' % (value['geo1+'],valenergy1,value['geo1/'],value['geo2+'],valenergy2,value['geo2/']))
                dertaenergy[num] ={}
                dertaenergy[num]['weight'] = value['weight']
                dertaenergy[num]['derta'] = math.fabs(valenergy - value['value'])
                dertaenergy[num]['QM-value'] = value['value']
                formenergy[num] = str(value['weight']).ljust(10) + formula.ljust(60) +\
                                  str('%.4f'% value['value']).rjust(10) +\
                                  str('%.4f'% valenergy).rjust(10) +\
                                  str('%.4f'% dertaenergy[num]['derta']).rjust(10)

            if numg == 3:
                formula = (value['geo1+'].ljust(2) + value['geo1'] + '/' + str(value['geo1/']).ljust(5)) +\
                          (value['geo2+'].ljust(2) + value['geo2'] + '/' + str(value['geo2/']).ljust(5)) +\
                          (value['geo3+'].ljust(2) + value['geo3'] + '/' + str(value['geo3/']).ljust(5))

                valenergy1 = fnomini[value['geo1']]['energy']
                valenergy2 = fnomini[value['geo2']]['energy']
                valenergy3 = fnomini[value['geo3']]['energy']
                valenergy = eval('%s%f/%f %s%f/%f %s%f/%f' % (value['geo1+'],valenergy1,value['geo1/'],\
                                                              value['geo2+'],valenergy2,value['geo2/'],\
                                                              value['geo3+'],valenergy3,value['geo3/']))
                dertaenergy[num] ={}
                dertaenergy[num]['weight'] = value['weight']
                dertaenergy[num]['derta'] = math.fabs(valenergy - value['value'])
                dertaenergy[num]['QM-value'] = value['value']
                formenergy[num] = str(value['weight']).ljust(10) + formula.ljust(60) +\
                                  str('%.4f'% value['value']).rjust(10) +\
                                  str('%.4f'% valenergy).rjust(10) +\
                                  str('%.4f'% dertaenergy[num]['derta']).rjust(10)

            if numg == 4:
                formula = (value['geo1+'].ljust(2) + value['geo1'] + '/' + str(value['geo1/']).ljust(5)) +\
                          (value['geo2+'].ljust(2) + value['geo2'] + '/' + str(value['geo2/']).ljust(5)) +\
                          (value['geo3+'].ljust(2) + value['geo3'] + '/' + str(value['geo3/']).ljust(5)) +\
                          (value['geo4+'].ljust(2) + value['geo4'] + '/' + str(value['geo4/']).ljust(5))

                valenergy1 = fnomini[value['geo1']]['energy']
                valenergy2 = fnomini[value['geo2']]['energy']
                valenergy3 = fnomini[value['geo3']]['energy']
                valenergy4 = fnomini[value['geo4']]['energy']
                valenergy = eval('%s%f/%f %s%f/%f %s%f/%f %s%f/%f' % (value['geo1+'],valenergy1,value['geo1/'],\
                                                                      value['geo2+'],valenergy2,value['geo2/'],\
                                                                      value['geo3+'],valenergy3,value['geo3/'],\
                                                                      value['geo4+'],valenergy4,value['geo4/']))
                dertaenergy[num] ={}
                dertaenergy[num]['weight'] = value['weight']
                dertaenergy[num]['derta'] = math.fabs(valenergy - value['value'])
                dertaenergy[num]['QM-value'] = value['value']
                formenergy[num] = str(value['weight']).ljust(10) + formula.ljust(60) +\
                                  str('%.4f'% value['value']).rjust(10) +\
                                  str('%.4f'% valenergy).rjust(10) +\
                                  str('%.4f'% dertaenergy[num]['derta']).rjust(10)

            if numg == 5:
                formula = (value['geo1+'].ljust(2) + value['geo1'] + '/' + str(value['geo1/']).ljust(5)) +\
                          (value['geo2+'].ljust(2) + value['geo2'] + '/' + str(value['geo2/']).ljust(5)) +\
                          (value['geo3+'].ljust(2) + value['geo3'] + '/' + str(value['geo3/']).ljust(5)) +\
                          (value['geo4+'].ljust(2) + value['geo4'] + '/' + str(value['geo4/']).ljust(5)) +\
                          (value['geo5+'].ljust(2) + value['geo5'] + '/' + str(value['geo5/']).ljust(5))

                valenergy1 = fnomini[value['geo1']]['energy']
                valenergy2 = fnomini[value['geo2']]['energy']
                valenergy3 = fnomini[value['geo3']]['energy']
                valenergy4 = fnomini[value['geo4']]['energy']
                valenergy5 = fnomini[value['geo5']]['energy']
                valenergy = eval('%s%f/%f %s%f/%f %s%f/%f %s%f/%f %s%f/%f' % (value['geo1+'],valenergy1,value['geo1/'],\
                                                                              value['geo2+'],valenergy2,value['geo2/'],\
                                                                              value['geo3+'],valenergy3,value['geo3/'],\
                                                                              value['geo4+'],valenergy4,value['geo4/'],\
                                                                              value['geo5+'],valenergy5,value['geo5/']))
                dertaenergy[num] ={}
                dertaenergy[num]['weight'] = value['weight']
                dertaenergy[num]['derta'] = math.fabs(valenergy - value['value'])
                dertaenergy[num]['QM-value'] = value['value']
                formenergy[num] = str(value['weight']).ljust(10) + formula.ljust(60) +\
                                  str('%.4f'% value['value']).rjust(10) +\
                                  str('%.4f'% valenergy).rjust(10) +\
                                  str('%.4f'% dertaenergy[num]['derta']).rjust(10)
    #
    # reaction_energy
    dertareacenergy = {}
    formreacenergy = {}
    if 'REACTION_ENERGY' in tset:

        for key,value in tset['REACTION_ENERGY'].items():
            num = key
            numg = value['numgeo']
            formula = []
            valenergy  = 0
            valenergy1 = 0
            valenergy2 = 0
            valenergy3 = 0
            valenergy4 = 0
            valenergy5 = 0
            if numg == 1:
                formula = value['geo1+'].ljust(2) + value['geo1'] + '/' + str(value['geo1/']).ljust(5)
                valenergy1 = fmini[value['geo1']]['energy']
                valenergy = valenergy1
                dertareacenergy[num] ={}
                dertareacenergy[num]['weight'] = value['weight']
                dertareacenergy[num]['derta'] = math.fabs(valenergy - value['value'])
                formreacenergy[num] = str(value['weight']).ljust(10) + formula.ljust(60) +\
                                  str('%.4f'% value['value']).rjust(10) +\
                                  str('%.4f'% valenergy).rjust(10) +\
                                  str('%.4f'% dertareacenergy[num]['derta']).rjust(10)

            if numg == 2:
                formula = (value['geo1+'].ljust(2) + value['geo1'] + '/' + str(value['geo1/']).ljust(5)) +\
                          (value['geo2+'].ljust(2) + value['geo2'] + '/' + str(value['geo2/']).ljust(5))

                valenergy1 = fmini[value['geo1']]['energy']
                valenergy2 = fmini[value['geo2']]['energy']
                valenergy = eval('%s%f/%f %s%f/%f' % (value['geo1+'],valenergy1,value['geo1/'],value['geo2+'],valenergy2,value['geo2/']))
                dertareacenergy[num] ={}
                dertareacenergy[num]['weight'] = value['weight']
                dertareacenergy[num]['derta'] = math.fabs(valenergy - value['value'])
                formreacenergy[num] = str(value['weight']).ljust(10) + formula.ljust(60) +\
                                  str('%.4f'% value['value']).rjust(10) +\
                                  str('%.4f'% valenergy).rjust(10) +\
                                  str('%.4f'% dertareacenergy[num]['derta']).rjust(10)

            if numg == 3:
                formula = (value['geo1+'].ljust(2) + value['geo1'] + '/' + str(value['geo1/']).ljust(5)) +\
                          (value['geo2+'].ljust(2) + value['geo2'] + '/' + str(value['geo2/']).ljust(5)) +\
                          (value['geo3+'].ljust(2) + value['geo3'] + '/' + str(value['geo3/']).ljust(5))

                valenergy1 = fmini[value['geo1']]['energy']
                valenergy2 = fmini[value['geo2']]['energy']
                valenergy3 = fmini[value['geo3']]['energy']
                valenergy = eval('%s%f/%f %s%f/%f %s%f/%f' % (value['geo1+'],valenergy1,value['geo1/'],\
                                                              value['geo2+'],valenergy2,value['geo2/'],\
                                                              value['geo3+'],valenergy3,value['geo3/']))
                dertareacenergy[num] ={}
                dertareacenergy[num]['weight'] = value['weight']
                dertareacenergy[num]['derta'] = math.fabs(valenergy - value['value'])
                formreacenergy[num] = str(value['weight']).ljust(10) + formula.ljust(60) +\
                                  str('%.4f'% value['value']).rjust(10) +\
                                  str('%.4f'% valenergy).rjust(10) +\
                                  str('%.4f'% dertareacenergy[num]['derta']).rjust(10)

            if numg == 4:
                formula = (value['geo1+'].ljust(2) + value['geo1'] + '/' + str(value['geo1/']).ljust(5)) +\
                          (value['geo2+'].ljust(2) + value['geo2'] + '/' + str(value['geo2/']).ljust(5)) +\
                          (value['geo3+'].ljust(2) + value['geo3'] + '/' + str(value['geo3/']).ljust(5)) +\
                          (value['geo4+'].ljust(2) + value['geo4'] + '/' + str(value['geo4/']).ljust(5))

                valenergy1 = fmini[value['geo1']]['energy']
                valenergy2 = fmini[value['geo2']]['energy']
                valenergy3 = fmini[value['geo3']]['energy']
                valenergy4 = fmini[value['geo4']]['energy']
                valenergy = eval('%s%f/%f %s%f/%f %s%f/%f %s%f/%f' % (value['geo1+'],valenergy1,value['geo1/'],\
                                                                      value['geo2+'],valenergy2,value['geo2/'],\
                                                                      value['geo3+'],valenergy3,value['geo3/'],\
                                                                      value['geo4+'],valenergy4,value['geo4/']))
                dertareacenergy[num] ={}
                dertareacenergy[num]['weight'] = value['weight']
                dertareacenergy[num]['derta'] = math.fabs(valenergy - value['value'])
                formreacenergy[num] = str(value['weight']).ljust(10) + formula.ljust(60) +\
                                  str('%.4f'% value['value']).rjust(10) +\
                                  str('%.4f'% valenergy).rjust(10) +\
                                  str('%.4f'% dertareacenergy[num]['derta']).rjust(10)

            if numg == 5:
                formula = (value['geo1+'].ljust(2) + value['geo1'] + '/' + str(value['geo1/']).ljust(5)) +\
                          (value['geo2+'].ljust(2) + value['geo2'] + '/' + str(value['geo2/']).ljust(5)) +\
                          (value['geo3+'].ljust(2) + value['geo3'] + '/' + str(value['geo3/']).ljust(5)) +\
                          (value['geo4+'].ljust(2) + value['geo4'] + '/' + str(value['geo4/']).ljust(5)) +\
                          (value['geo5+'].ljust(2) + value['geo5'] + '/' + str(value['geo5/']).ljust(5))

                valenergy1 = fmini[value['geo1']]['energy']
                valenergy2 = fmini[value['geo2']]['energy']
                valenergy3 = fmini[value['geo3']]['energy']
                valenergy4 = fmini[value['geo4']]['energy']
                valenergy5 = fmini[value['geo5']]['energy']
                valenergy = eval('%s%f/%f %s%f/%f %s%f/%f %s%f/%f %s%f/%f' % (value['geo1+'],valenergy1,value['geo1/'],\
                                                                              value['geo2+'],valenergy2,value['geo2/'],\
                                                                              value['geo3+'],valenergy3,value['geo3/'],\
                                                                              value['geo4+'],valenergy4,value['geo4/'],\
                                                                              value['geo5+'],valenergy5,value['geo5/']))
                dertareacenergy[num] ={}
                dertareacenergy[num]['weight'] = value['weight']
                dertareacenergy[num]['derta'] = math.fabs(valenergy - value['value'])
                formreacenergy[num] = str(value['weight']).ljust(10) + formula.ljust(60) +\
                                  str('%.4f'% value['value']).rjust(10) +\
                                  str('%.4f'% valenergy).rjust(10) +\
                                  str('%.4f'% dertareacenergy[num]['derta']).rjust(10)


    return(dertacharge,dertageometry,dertaenergy,dertareacenergy,dertacellpam,dertaheat,formcharge,formgeometry,formenergy,formreacenergy,formcellpam,formheat)
