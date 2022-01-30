# -*- coding:utf-8 -*-
import os,sys
from lammps import lammps

def ybs_lmpnomini(ffi,nom,restraint,geonatom,datatype,cell,lg,m=0,outdir='outdata'):
    #  can use mpi
    chirld_geo = nom
    ffield = ffi
    mi = m
    fnomini = {}
    for key in chirld_geo:

        geo = key
        dataname = "data." + key
        lmp = lammps()
        lmp.command("units real")
        lmp.command("dimension       3")
        #
        if geo in cell:
            lmp.command("boundary        p p p")
        else:
            lmp.command("boundary        f f f")
        lmp.command("atom_style charge")
        if (geo in restraint) and mi == 1:
            lmp.command("atom_modify map array")
        lmp.command("read_data ./data/%s" % dataname)
        if 1 == lg:
            lmp.command("pair_style reax/c NULL safezone 10 mincap 2000 lgvdw yes")
        else:
            lmp.command("pair_style reax/c NULL safezone 10 mincap 2000")
        lmp.command("pair_coeff * * ./%s/%s %s" % (outdir, ffield, datatype[dataname]))
        # if geo in cell:
        lmp.command("neighbor	2.5 bin")
        lmp.command("neigh_modify	every 1 delay 0 check yes one 2000 page 100000")
        lmp.command("fix             1 all qeq/reax 1 0.0 10.0 1.0e-6 reax/c")
        lmp.command("timestep 0.1")
        if mi == 1:
            if geo in restraint:
                lmp.command("fix hold all %s" % restraint[geo])
                lmp.command("fix_modify hold energy yes")
                lmp.command("min_style fire")
            else:
                lmp.command("min_style cg")
            if geo in cell:
                pass
                #lmp.command("fix 2 all box/relax iso 0.0 vmax 0.001")
            # lmp.command("min_modify line forcezero")
            lmp.command("minimize 0.0 1.0e-4 1000 10000")
        lmp.command("run 0")
        lmp.command("compute reax all pair reax/c")
        qq = lmp.gather_atoms('q', 1, 1)
        q = [0 for i in range(geonatom[geo])]
        q = [qq[i] for i in range(geonatom[geo])]
        xx = lmp.gather_atoms('x', 1, 3)
        x1 = [0 for i in range(3 * geonatom[geo])]
        x1 = [xx[i] for i in range(3 * geonatom[geo])]
        x = ([x1[i:i + 3] for i in range(0, len(x1), 3)])
        eng = lmp.extract_compute("reax", 0, 0)
        fnomini[geo] = {}
        fnomini[geo]['q'] = q
        fnomini[geo]['x'] = x
        fnomini[geo]['energy'] = eng
        lmp.close()
    return fnomini
##########################################################################################
def ybs_lmpmini(ffi, nom,restraint,geonatom,datatype,cell,lg,m=0,outdir='outdata'):
    #  can use mpi
    chirld_geo = nom
    ffield = ffi
    mi = m
    fmini = {}
    if mi ==3:
        for key in chirld_geo:
            geo = key
            dataname = "data." + key
            try:
                lmp = lammps()
                lmp.command("units real")
                lmp.command("dimension       3")
                #
                if geo in cell:
                    lmp.command("boundary        p p p")
                else:
                    lmp.command("boundary        f f f")
                lmp.command("atom_style charge")
                #if geo in restraint:
                    #lmp.command("atom_modify map array")
                lmp.command("read_data ./data/%s" % dataname)
                if 1 == lg:
                    lmp.command("pair_style reax/c NULL safezone 10 mincap 2000 lgvdw yes")
                else:
                    lmp.command("pair_style reax/c NULL safezone 10 mincap 2000")
                lmp.command("pair_coeff * * ./%s/%s %s" % (outdir,ffield, datatype[dataname]))
                #if geo in cell:
                lmp.command("neighbor	2.5 bin")
                lmp.command("neigh_modify	every 1 delay 0 check yes one 2000 page 100000")
                lmp.command("fix             1 all qeq/reax 1 0.0 10.0 1.0e-6 reax/c")
                lmp.command("timestep 0.1")
                #if geo in restraint:
                    #lmp.command("fix hold all %s" % restraint[geo])
                    #lmp.command("fix_modify hold energy yes")
                    #lmp.command("min_style fire")
                #else:
                    #lmp.command("min_style cg")
                #if geo in cell:
                    #lmp.command("fix 2 all box/relax iso 0.0 vmax 0.001")
                #lmp.command("min_modify line forcezero")
                #lmp.command("minimize 0.0 1.0e-4 1000 10000")
                lmp.command("thermo_style custom step cella cellb cellc")
                lmp.command("run 0")
                lmp.command("compute reax all pair reax/c")
                qq = lmp.gather_atoms('q', 1, 1)
                q = [0 for i in range(geonatom[geo])]
                q = [qq[i] for i in range(geonatom[geo])]
                xx = lmp.gather_atoms('x', 1, 3)
                x1 = [0 for i in range(3 * geonatom[geo])]
                x1 = [xx[i] for i in range(3 * geonatom[geo])]
                x = ([x1[i:i + 3] for i in range(0, len(x1), 3)])
                eng = lmp.extract_compute("reax", 0, 0)
                fmini[geo] = {}
                fmini[geo]['q'] = q
                fmini[geo]['x'] = x
                fmini[geo]['energy'] = eng
                if geo in cell:
                    a = lmp.get_thermo('cella')
                    b = lmp.get_thermo('cellb')
                    c = lmp.get_thermo('cellc')
                    density = lmp.get_thermo('density')
                    press = lmp.get_thermo('press')
                    fmini[geo]['a'] = a
                    fmini[geo]['b'] = b
                    fmini[geo]['c'] = c
                    fmini[geo]['r'] = density
                    fmini[geo]['p'] = press
                lmp.close()
            except:
                print(geo)
                sys.exit()
    else:
        for key in chirld_geo:
            geo = key
            dataname = "data." + key
            lmp = lammps()
            lmp.command("units real")
            lmp.command("dimension       3")
            if geo in cell:
                lmp.command("boundary        p p p")
            else:
                lmp.command("boundary        p p p")
            lmp.command("atom_style charge")
            if geo in restraint:
                lmp.command("atom_modify map array")
            lmp.command("read_data ./data/%s" % dataname)
            if 1 == lg:
                lmp.command("pair_style reax/c NULL lgvdw yes")
            else:
                lmp.command("pair_style reax/c NULL")
            lmp.command("pair_coeff * * ./%s/%s %s" % (outdir,ffield,datatype[dataname]))
            #if geo in cell:
            lmp.command("neighbor	2.5 bin")
            lmp.command("neigh_modify	every 1 delay 0 check yes one 2000 page 100000")
            lmp.command("fix             1 all qeq/reax 1 0.0 10.0 1.0e-6 reax/c")
            lmp.command("timestep 0.2")
            if geo in restraint:
                lmp.command("fix hold all %s" % restraint[geo])
                lmp.command("fix_modify hold energy yes")
                lmp.command("min_style fire")
            else:
                lmp.command("min_style cg")
            if geo in cell:
                lmp.command("fix 2 all box/relax iso 0.0 vmax 0.001")
            #lmp.command("min_modify line forcezero")
            lmp.command("thermo_style custom step density press cella cellb cellc")
            lmp.command("minimize 1.0e-9 1.0e-6 200 1000")
            #lmp.command("run 0")
            #lmp.command("fix             10knvt all nvt temp 10.0 10.0 1")
            #lmp.command("run             5000")
            #lmp.command("unfix           10knvt")
            if 1 == 2:
                print(geo)
                lmp.command("fix             300knvt all nvt temp 300.0 300.0 10")
                lmp.command("run             2000")
                lmp.command("unfix           300knvt")
                lmp.command("fix             0knvt all nvt temp 300.0 100.0 1")
                lmp.command("run             500")
                lmp.command("unfix           0knvt")
                lmp.command("minimize 0.0 1.0e-6 1000 10000")
            lmp.command("compute reax all pair reax/c")
            qq = lmp.gather_atoms('q',1,1)
            q = [0 for i in range(geonatom[geo])]
            q = [qq[i] for i in range(geonatom[geo])]
            xx = lmp.gather_atoms('x',1,3)
            x1 = [0 for i in range(3*geonatom[geo])]
            x1 = [xx[i] for i in range(3*geonatom[geo])]
            x  = ([x1[i:i + 3] for i in range(0, len(x1), 3)])
            eng = lmp.extract_compute("reax",0,0)
            fmini[geo] = {}
            fmini[geo]['q'] = q
            fmini[geo]['x'] = x
            fmini[geo]['energy'] = eng
            if geo in cell:
                a = lmp.get_thermo('cella')
                b = lmp.get_thermo('cellb')
                c = lmp.get_thermo('cellc')
                density = lmp.get_thermo('density')
                press = lmp.get_thermo('press')
                fmini[geo]['a'] = a
                fmini[geo]['b'] = b
                fmini[geo]['c'] = c
                fmini[geo]['r'] = density
                fmini[geo]['p'] = press
            lmp.close()
    return fmini
