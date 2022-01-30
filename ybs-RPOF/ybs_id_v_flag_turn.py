import os
import math, re


def id_v_flag_turn(pam,new_pam_id_v_flag):
    from ybs_ffieldpam import ffpam
    # start pam flag
    x_flag_value_str_list = ""
    for num in pam:
        pam_a = pam[num]['value'][1]
        pam_b = pam[num]['value'][2]
        id1 = int(pam[num]['id'][1])
        id2 = int(pam[num]['id'][2])
        id3 = int(pam[num]['id'][3])
        pam_v = ffpam("ffield.reax.start", [id1, id2, id3])
        pam_flag = str(id1) + '-' + str(id2) + '-' + str(id3)
        x_flag_value_str = pam_flag + "_" + str(pam_v) + "_"
        x_flag_value_str_list += x_flag_value_str
    # orignal
    flag_list1 = x_flag_value_str_list.split("_")
    flag_list2 = [flag_list1[i:i+2] for i in range(0,len(flag_list1),2)]
    flag_dir = {}
    for i in flag_list2:
        if i[0] != "":
            flag_dir[i[0]] = i[1]

    # new
    flag_list3 = new_pam_id_v_flag.split("_")
    flag_list4 = [flag_list3[i:i+2] for i in range(0,len(flag_list3),2)]
    flag_dir2 = {}
    for i in flag_list4:
        if i[0] != "":
            flag_dir2[i[0]] = i[1]

    # change pam_v
    for k,v in flag_dir2.items():
        pam_flag = k
        pam_v = v
        flag_dir[pam_flag] = pam_v

    x_flag_value_str_list_new = ""
    for k,v in flag_dir.items():
        pam_flag = k
        pam_v = v
        #x_flag_value_str = pam_flag + "_" + str(pam_v) + "_"
        x_flag_value_str = str(pam_v)
        x_flag_value_str_list_new += x_flag_value_str
    return x_flag_value_str_list,x_flag_value_str_list_new

if __name__ == "__main__":
    with open("params") as r:
        lines = r.readlines()
    pam = {}
    numpam = 0
    for i in lines:
        perpam = i.split()
        if len(perpam) > 0 and i.find("#") != 0:
            numpam += 1
            pam[numpam] = {}
            pam[numpam]['id'] = {}
            pam[numpam]['value'] = {}
            pam[numpam]['id'][1] = int(perpam[0])
            pam[numpam]['id'][2] = int(perpam[1])
            pam[numpam]['id'][3] = int(perpam[2])
            pam[numpam]['value'][1] = float(perpam[-2])
            pam[numpam]['value'][2] = float(perpam[-1])
            pam[numpam]['step'] = float(perpam[-3])
    new_pam_id_v_flag = "4-6-1_0.095"
    a = id_v_flag_turn(pam,new_pam_id_v_flag)
