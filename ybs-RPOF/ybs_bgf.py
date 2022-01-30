def bgfdata(inputgeo):
    import os, sys
    import  string
    from collections import Counter
    import numpy as np
    from math import sqrt,fabs,cos,sin,pi
    numgeos = 0
    dir_path = os.getcwd()
    bgf = open(dir_path + os.sep + inputgeo,'r')
    bgflines=bgf.readlines()
    datanametype ={}
    geos={}
    types={}
    atoms=0
    numAtoms={}
    restraintname ={}
    cell = {}
    yuanzi = {}
    geontype = {}
    for line in bgflines:
        if line.find('DESCRP')==0:
            name=line.split()[1]
            geos[name]={}
            types[name]=[]
            yuanzi[name] = []
            numgeos += 1
            #print(line,name)
        elif line.find('END')==0:
            numAtoms[name]=atoms
            atoms=0
        elif line.find('HETATM')==0:
            atoms+=1
            tokens=line.split()
            id=int(tokens[1])
            #print name,id,type,x,y,z
            geos[name][id]={}
            geos[name][id]['type']=tokens[6]
            geos[name][id]['x']=float(tokens[3])
            geos[name][id]['y']=float(tokens[4])
            geos[name][id]['z']=float(tokens[5])
            yuanzi[name].append(tokens[6])
            if tokens[6] not in types[name]:
                types[name].append(tokens[6])
        elif line.find('BOND')==0:
            resline = line.split()
            restraintname[name] = "restrain bond %d %d 8000.0 8000.0 %f" % (int(resline[2]),int(resline[3]),float(resline[4]))
        elif line.find('ANGLE')==0:
            resline = line.split()
            if len(resline) == 9: # 
                restraintname[name] = "restrain angle %d %d %d 7000.0 7000.0 %f" % (int(resline[2]),int(resline[3]),int(resline[4]),float(resline[5]))
            if len(resline) == 10: # 
                if float(resline[6]) > 180:
                    restrain_dihedral_value = float(resline[6]) - 360.0
                else:
                    restrain_dihedral_value = float(resline[6])
                restraintname[name] = "restrain dihedral %d %d %d %d 6000.0 6000.0 %f" % (
                int(resline[2]), int(resline[3]), int(resline[4]), int(resline[5]), (restrain_dihedral_value + 180))
        elif line.find('DIHEDRAL')==0:
            resline = line.split()
            if float(resline[6]) > 180:
                restrain_dihedral_value = float(resline[6]) - 360.0
            else:
                restrain_dihedral_value = float(resline[6])
            restraintname[name] = "restrain dihedral %d %d %d %d 6000.0 6000.0 %f" % (
                int(resline[2]), int(resline[3]), int(resline[4]), int(resline[5]), (restrain_dihedral_value + 180))
        elif line.find('CRYSTX')==0:
            cell[name] = {}
            cell[name]['type'] = 'cell'
            cell[name]['abc'] = line.split()[1:]
        else:
            continue
    for key,value in yuanzi.items():
        geo =key
        geontype[geo]=[]
        geontype[geo]=Counter(yuanzi[geo])

    for key,value in cell.items():
        geo = key
        A = float(value['abc'][0])
        B = float(value['abc'][1])
        C = float(value['abc'][2])
        alpha = float(value['abc'][3])
        beta  = float(value['abc'][4])
        gamma = float(value['abc'][5])
        eps = 2 * np.spacing(90.0, dtype=np.float64)

        # alpha
        if abs(abs(alpha) - 90) < eps:
            cos_alpha = 0.0
        else:
            cos_alpha = cos(alpha * pi / 180.0)
        # beta
        if abs(abs(beta) - 90) < eps:
            cos_beta = 0.0
        else:
            cos_beta = cos(beta * pi / 180.0)
        # gamma
        if abs(gamma - 90) < eps:
            cos_gamma = 0.0
            sin_gamma = 1.0
        elif abs(gamma + 90) < eps:
            cos_gamma = 0.0
            sin_gamma = -1.0
        else:
            sin_gamma = sin(gamma * pi / 180.0)
            cos_gamma = cos(gamma * pi / 180.0)
        # Build the cell vectors
        cx = cos_beta
        cy = (cos_alpha - cos_beta * cos_gamma) / sin_gamma
        cz = sqrt(1. - cx * cx - cy * cy)
        va = A * np.array([1, 0, 0])
        vb = B * np.array([cos_gamma, sin_gamma, 0])
        vc = C * np.array([cx, cy, cz])
        # Convert to the Cartesian x,y,z-system
        abc = np.vstack((va, vb, vc))
        # T = np.vstack((X, Y, Z))
        # cell = dot(abc, T)

        a = va[0]
        b = vb[1]
        c = vc[2]
        xy = vb[0]
        xz = vc[0]
        yz = vc[1]
        cell[geo]['x'] = '%.6f'% a
        cell[geo]['y'] = '%.6f'% b
        cell[geo]['z'] = '%.6f'% c
        cell[geo]['xy'] = '%.6f'% xy
        cell[geo]['xz'] = '%.6f'% xz
        cell[geo]['yz'] = '%.6f'% yz
    if not os.path.exists(os.curdir+'/'+"data"):
        os.mkdir(os.curdir+os.sep+"data")
    else:
        pass
    for key,value in geos.items():
        geoname = key
        geoatoms = value
        atomtypes = types[geoname]
        atomtypes.sort()
        numggeoatom = numAtoms[geoname]
        numatomstr = str(numggeoatom)
        typenum = len(types[geoname])
        typenumstr = str(typenum)
        dataname = "data."+ geoname
        datanametype[dataname] = ' '.join(atomtypes)
        datafilepath = dir_path+os.sep+"data"+os.sep+dataname
        out = open(datafilepath,'w')
        out.write("")
        out.close()
        out = open(datafilepath,'a')
        out.write(geoname + '\n' + numatomstr + " atoms" + '\n')
        out.write(typenumstr + " atom types" + '\n' + '\n')
        lo = '%.6f' % (geos[geoname][1]['x'] - 100)
        hi = '%.6f' % (geos[geoname][1]['x'] + 100)
        if geoname not in  cell:
            out.write(str(lo).rjust(12) + str(hi).rjust(12) + " xlo" + " xhi" + '\n')
            out.write(str(lo).rjust(12) + str(hi).rjust(12) + " ylo" + " yhi" + '\n')
            out.write(str(lo).rjust(12) + str(hi).rjust(12) + " zlo" + " zhi" + '\n' + '\n')
        else:
            out.write(str('%.6f'% 0.0).rjust(12) + str(cell[geoname]['x']).rjust(12) + " xlo" + " xhi" + '\n')
            out.write(str('%.6f'% 0.0).rjust(12) + str(cell[geoname]['y']).rjust(12) + " ylo" + " yhi" + '\n')
            out.write(str('%.6f'% 0.0).rjust(12) + str(cell[geoname]['z']).rjust(12) + " zlo" + " zhi" + '\n')
            if cell[geo]['xy'] == cell[geo]['xz'] == cell[geo]['yz'] == 0.:
                pass
            else:
                out.write(str(cell[geoname]['xy']).rjust(12) + str(cell[geoname]['xz']).rjust(12) + str(cell[geoname]['yz']).rjust(12) +" xy" +" xz" + " yz" + '\n' + '\n')
        out.write("Masses" + '\n' + '\n')
        for i in atomtypes:
            typename = i
            typeid  = atomtypes.index(typename) + 1
            if typename == 'C':
                out.write(str(typeid) + ' ' + '12.011000' + ' # C' + '\n')
            elif typename == 'H':
                out.write(str(typeid) + ' ' + '1.008000' + '  # H' + '\n')
            elif typename == 'N':
                out.write(str(typeid) + ' ' + '14.007000' + ' # N' + '\n')
            elif typename == 'O':
                out.write(str(typeid) + ' ' + '15.999000' + ' # O' + '\n')
            elif typename == 'Al':
                out.write(str(typeid) + ' ' + '26.982000' + ' # Al' + '\n')
            elif typename == 'Cl':
                out.write(str(typeid) + ' ' + '35.450000' + ' # Cl' + '\n')
            elif typename == 'Si':
                out.write(str(typeid) + ' ' + '28.085500' + ' # Si' + '\n')
            else:
                print(i)
        out.write('\n' + "Atoms" + '\n' + '\n')
        for key1,value1 in (geoatoms.items()):
            atomid = key1
            atomdata = value1
            typeid = atomtypes.index(atomdata['type']) + 1
            x =  '%.6f' % atomdata['x']
            y =  '%.6f' % atomdata['y']
            z =  '%.6f' % atomdata['z']
            out.write(str(atomid) + '  ' + str(typeid) + '  0.000000 ' + str(x).rjust(15) + str(y).rjust(15) +str(z).rjust(15) + ' # '+ atomdata['type'] + '\n')
        out.close()
    return(numAtoms,numgeos,datanametype,restraintname,cell,geontype)

