# -*- coding:utf-8 -*-
import os,sys

def reassign_geo(mpi_nc,mnp,mini,nomini,cell,md=1):
    mpi_m = mnp
    mpi_n = mpi_nc - mpi_m
    if md ==1:
        mpi_n = mpi_nc
        mpi_m = 0
        if len(mini) != 0:
            if mpi_nc > len(mini):
                mpi_n = mpi_nc - len(mini)
                mpi_m = len(mini)
                #print(mpi_nc,mpi_n,mpi_m)
            else:
                sys.exit()
    #
    if len(nomini) != 0:
        chirld_nomini_cell = []
        cell_nomini_list = []
        for geo_name in nomini:
            if geo_name in cell and "NTO-f" not in geo_name:
                cell_nomini_list.append(geo_name)
                #nomini.remove(geo_name)
        for geo_name in cell_nomini_list:
            #print(geo_name)
            if geo_name in nomini:
                nomini.remove(geo_name)
                #print("yes")
        #
        if len(nomini) % mpi_n == 0:
            chirldnonum = (len(nomini) // mpi_n)
            chirld_nomini = [nomini[i:i + chirldnonum] for i in range(0, len(nomini), chirldnonum)]
        else:
            #
            chirldnonum = (len(nomini) // mpi_n)
            chirld_nomini = [nomini[i:i + chirldnonum] for i in range(0, chirldnonum * mpi_n, chirldnonum)]
            yu_nomini = nomini[chirldnonum * mpi_n:]
            for i in range(len(yu_nomini)):
                chirld_nomini[i].append(yu_nomini[i])
        if 0 < len(cell_nomini_list) <= mpi_n:
            for i in range(len(cell_nomini_list)):
                chirld_nomini[i] += [cell_nomini_list[i]]
        elif mpi_n < len(cell_nomini_list):
            ####
            if len(cell_nomini_list) != 0:
                if len(cell_nomini_list) % mpi_n == 0:
                    chirldnonum_cell = (len(cell_nomini_list) // mpi_n)
                    chirld_nomini_cell = [cell_nomini_list[i:i + chirldnonum_cell] for i in range(0, len(cell_nomini_list), chirldnonum_cell)]
                # cell
                chirldnonum_cell = (len(cell_nomini_list) // mpi_n)
                chirld_nomini_cell = [cell_nomini_list[i:i + chirldnonum_cell] for i in range(0, chirldnonum_cell * mpi_n, chirldnonum_cell)]
                yu_nomini_cell = cell_nomini_list[chirldnonum_cell * mpi_n:]
                for i in range(len(yu_nomini_cell)):
                    chirld_nomini_cell[i].append(yu_nomini_cell[i])
            for i in range(mpi_n):
                chirld_nomini[i] += chirld_nomini_cell[i]
    # mini
    if len(mini) != 0:
        #
        if len(mini) % mpi_m == 0:
            chirldnum = (len(mini) // mpi_m)
            chirld_mini = [mini[i:i + chirldnum] for i in range(0, len(mini), chirldnum)]
        else:
            #
            chirldnum = (len(mini) // mpi_m)
            chirld_mini = [mini[i:i + chirldnum] for i in range(0, chirldnum * mpi_m, chirldnum)]
            yu_mini = mini[chirldnum * mpi_m:]
            for i in range(len(yu_mini)):
                chirld_mini[i].append(yu_mini[i])
    return chirld_mini,chirld_nomini,chirld_nomini_cell,mpi_n,mpi_m
