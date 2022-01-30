import math
#bond
def bond(at1,at2):
    atom1 = at1
    atom2 = at2
    bonds = 0
    for i in range(3):
        bonds += math.pow((atom1[i]-atom2[i]),2)
    bonds = math.sqrt(bonds)
    return(bonds)
def angle(at1,at2,at3):
    #  angle
    atom1 = at1
    atom2 = at2
    atom3 = at3
    dx1 = atom1[0] - atom2[0]
    dy1 = atom1[1] - atom2[1]
    dz1 = atom1[2] - atom2[2]
    rsq1 = dx1 * dx1 + dy1 * dy1 + dz1 * dz1
    r1 = math.sqrt(rsq1)
    #
    dx2 = atom3[0] - atom2[0]
    dy2 = atom3[1] - atom2[1]
    dz2 = atom3[2] - atom2[2]
    rsq2 = dx2 * dx2 + dy2 * dy2 + dz2 * dz2
    r2 = math.sqrt(rsq2)
    #
    c = dx1 * dx2 + dy1 * dy2 + dz1 * dz2
    c /= r1 * r2
    if c > 1:
        c = 1.0
    elif c < -1.0:
        c= -1.0
    theta = math.acos(c)
    jiao = theta * 180 / math.pi
    return(jiao)
#
def torsion(at1,at2,at3,at4):
    atom1 = at1
    atom2 = at2
    atom3 = at3
    atom4 = at4
    ip00 = ip11 = ip22 = ip01 = ip02 = ip12 = ip01x2 = 0.0
    v=[[0 for i in range(3)] for i in range(3)]
    n = [0 for i in range(3)]
    for k in range(3):
        v[0][k] = (atom2[k]) - (atom1[k])
        v[1][k] = (atom3[k]) - (atom2[k])
        v[2][k] = (atom4[k]) - (atom3[k])
        ip00 += (v[0][k] * v[0][k])
        ip11 += (v[1][k] * v[1][k])
        ip22 += (v[2][k] * v[2][k])
        ip01 += (v[0][k] * v[1][k])
        ip02 += (v[0][k] * v[2][k])
        ip12 += (v[1][k] * v[2][k])
    i = 1
    j = 2
    n[0] = v[i][1] * v[j][2] - v[i][2] * v[j][1]
    n[1] = v[i][2] * v[j][0] - v[i][0] * v[j][2]
    n[2] = v[i][0] * v[j][1] - v[i][1] * v[j][0]
    for k in range(3):
        ip01x2 += (v[0][k] * n[k])
        #x=-v[0]*v[2]+(v[0]*v[1])(v[1]*v[2])= v[0]*(-v[2]+(v[1]*v[2])v[1])
        x = -ip02 / math.sqrt(ip00 * ip22) + (ip01 * ip12) / (ip11 * math.sqrt(ip00 * ip22))
        #y=v[0]*(v[1]xv[2])
        y = -ip01x2 / math.sqrt(ip00 * ip11 * ip22)
        #printf("%f %f %f %f %f %f %f \n",ip00,ip11,ip22,ip01,ip02,ip12,ip01x2)
        #printf("%f %f \n",x,y)
    costh = x / math.sqrt(x * x + y * y)
    if costh > 1:
        costh = 1.0
    elif costh < -1.0:
        costh = -1.0
    theta = math.acos(costh)
    if (y > 0):
        theta *= -1
    jiao2 = theta * 180 / math.pi
    return(jiao2)

