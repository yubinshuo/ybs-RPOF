import math
def Fx(x,f,t,w1=1.0,w2=1.0,w3=1.0,w4=1.0, w5=1.0, w6=1.0,fix={}):
    ffieldname = x
    tset = t
    Fx = 0
    fix_error = fix
    setout = {}
    sumcharge = 0
    sumgeometry = 0
    sumenergy = 0
    # reaction_energy
    sumreacenergy = 0
    sumcellpam = 0
    sumheat = 0
    dcharge,dgeometry,denergy,dreacenergy,dcellpam,dheat,formcha,formgeo,formene,formreacene,formcell,formheat = f
    #
    numcha = len(formcha)
    numgeo = len(formgeo)
    numene = len(formene)
    numreacene = len(formreacene)
    numcell = len(formcell)
    numheat = len(formheat)
    sumentries = numheat+numene+numcha+numcell+numgeo+numreacene
    sumsecw = (numcha > 0)*1 + (numgeo > 0)*1 +(numene > 0)*1 + (numcell > 0)*1 +(numheat > 0)*1 + (numreacene > 0)*1
    wseccha = 1
    wsecgeo = 1
    wsecene = 1
    wsecreacene = 1
    wseccell = 1
    wsecheat = 1
    if w1 ==1 and w2==1 and w3==1 and w4==1 and w5==1and w6==1:
        wseccha = numcha / sumentries
        wsecgeo = numgeo / sumentries
        wsecene = numene / sumentries
        wsecreacene = numreacene / sumentries
        wseccell = numcell / sumentries
        wsecheat = numheat / sumentries
    else:
        sumsecweight = w1 + w2 + w3 + w4 + w5 + w6
        wseccha = w1 / sumsecweight
        wsecgeo = w2 / sumsecweight
        wsecene = w3 / sumsecweight
        wsecreacene = w4 / sumsecweight
        wseccell = w5 / sumsecweight
        wsecheat = w6 / sumsecweight
    ww = 0
    if w1 ==2 and w2==2 and w3==2 and w4==2 and w5==2 and w6==2:
        ww = 1

    for key,value in formcha.items():
        num = key
        setout[num] = value
    for key,value in formgeo.items():
        num = len(formcha) + key
        setout[num] = value
    for key,value in formene.items():
        num = len(formcha) + len(formgeo) + key
        setout[num] = value
    # reaction_energy
    for key,value in formreacene.items():
        num = len(formcha) + len(formgeo) + len(formene) + key
        setout[num] = value

    for key,value in formcell.items():
        num = len(formcha) + len(formgeo) + len(formene) + len(formreacene)+ key
        setout[num] = value
    for key,value in formheat.items():
        num = len(formcha) + len(formgeo) + len(formene)+ len(formreacene) + len(formcell) + key
        setout[num] = value
    ss,s1,s2,s3,s4,s5,s6  = 0,0,0,0,0,0,0

    #charge
    if 'CHARGE' in tset:
        sumchawi = 0.0
        for key,value in dcharge.items():
            wi = value['weight']
            sumchawi += wi
        for key,value in dcharge.items():
            fc = value['derta']
            wi = value['weight']
            weight = (wi / sumchawi)*wseccha
            sigma = wi
            if ww == 1:
                sumcharge += ((fc / wi) * (fc / wi))
            elif ww == 2:
                sumcharge += ((fc / wi) * (fc / wi))
            else:
                sumcharge += ((fc/wi) * (fc/wi))*wseccha
            s1 += weight
            ss += weight
    if 'GEOMETRY' in tset:
        # geometry
        sumgeowi = 0.0
        for key, value in dgeometry.items():
            wi = value['weight']
            sumgeowi += wi
        for key,value in dgeometry.items():
            num = key
            fc = value['derta']
            wi = value['weight']
            if tset['GEOMETRY'][num]['type']=='bond':
                #
                sigma_bond_w1 = 0.1
                sigma_bond_w2 = 8
                sigma = sigma_bond_w1 * math.pow(sigma_bond_w1,(sigma_bond_w2 * fc)) + 0.0001
                #
            elif tset['GEOMETRY'][num]['type']=='angle' or tset['GEOMETRY'][num]['type']=='torsion':
                #
                sigma_ad_w1 = 1.5
                sigma_ad_w2 = 0.4
                sigma = sigma_ad_w1 * math.pow((0.1 * sigma_ad_w1),(sigma_ad_w2 * fc)) + 0.5
                #
            weight = (wi / sumgeowi)*wsecgeo
            if ww == 1:
                sumgeometry += ((fc / wi) * (fc / wi))
            elif ww == 2:
                sumgeometry += ((fc / wi) * (fc / wi))
            else:
                sumgeometry += ((fc / wi) * (fc / wi))*wsecgeo
            s2 += weight
            ss += weight
    # energy
    if 'ENERGY' in tset:
        sumenewi = 0.0
        for key, value in denergy.items():
            wi = value['weight']
            sumenewi += wi
        for key,value in denergy.items():

            fc = value['derta']
            wi = value['weight']
            qm_value = value['QM-value']
            entry_num = key + numcha + numgeo
            weight = (wi / sumenewi)*wsecene
            sigma_energy_w1 = 4.0
            sigma_energy_w2 = 0.5
            sigma_energy_w3 = 0.02
            sigma_energy_w = math.fabs(qm_value) + math.e
            sigma = -sigma_energy_w1 * math.pow( sigma_energy_w2, (sigma_energy_w3 * sigma_energy_w)) + sigma_energy_w1 + 0.1

            if entry_num in fix_error:
                sigma = 0.4 * sigma + 0.1
            if ww == 1:
                sumenergy += ((fc / wi) * (fc / wi))
            elif ww == 2:
                sumenergy += ((fc / wi) * (fc / wi))
            else:
                sumenergy += ((fc / wi) * (fc / wi))*wsecene
            s3 += weight
            ss += weight
    if 'REACTION_ENERGY' in tset:
        sumreacenewi = 0.0
        for key,value in dreacenergy.items():
            wi = value['weight']
            sumreacenewi += wi
        for key,value in dreacenergy.items():
            fc = value['derta']
            wi = value['weight']
            weight = (wi / sumreacenewi)*wsecreacene
            if ww == 1:
                sumreacenergy += ((fc / wi) * (fc / wi))
            elif ww == 2:
                sumreacenergy += ((fc / wi) * (fc / wi))
            else:
                sumreacenergy += ((fc / wi) * (fc / wi)) * wsecreacene
            s4 += weight
            ss += weight
    if 'CELL_PARAMETERS' in tset:
        sumcellwi = 0.0
        for key,value in dcellpam.items():
            wi = value['weight']
            sumcellwi += wi
        for key,value in dcellpam.items():
            fc = value['derta']
            wi = value['weight']
            weight = (wi /sumcellwi)*wseccell
            sigma = 0.1
            if ww == 1:
                sumcellpam += ((fc / wi) * (fc / wi))
            elif ww == 2:
                sumcellpam += ((fc / wi) * (fc / wi)) * wi
            else:
                sumcellpam += ((fc / wi) * (fc / wi))*wseccell
            s5 += weight
            ss += weight

    # heat
    if 'HEATFO' in tset:
        sumheatwi = 0.0
        for key,value in dheat.items():
            wi = value['weight']
            sumheatwi += wi
        for key,value in dheat.items():
            fc = value['derta']
            wi = value['weight']
            qm_value = value['QM-value']
            weight = (wi /sumheatwi)*wsecheat
            sigma_heat_w1 = 4.0
            sigma_heat_w2 = 0.5
            sigma_heat_w3 = 0.02
            sigma_heat_w = math.fabs(qm_value) + math.e
            sigma = -sigma_heat_w1 * math.pow(sigma_heat_w2, (sigma_heat_w3 * sigma_heat_w)) + sigma_heat_w1 + 0.1
            if ww == 1:
                sumheat += ((fc / wi) * (fc / wi))
            elif ww == 2:
                sumheat += ((fc / wi) * (fc / wi))
            else:
                sumheat += ((fc / wi) * (fc / wi))*wsecheat
            s6 += weight
            ss += weight
    Fx = sumcharge + sumgeometry + sumenergy + sumcellpam + sumheat + sumreacenergy
    return(Fx,sumcharge,sumgeometry,sumenergy,sumreacenergy,sumcellpam,sumheat,setout)

