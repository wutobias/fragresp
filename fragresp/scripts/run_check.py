import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pickle
from collections import OrderedDict

from fragresp.gaussian_tools import check_opt as _check_opt
from fragresp.gaussian_tools import check_esp as _check_esp
from fragresp.gaussian_tools import get_energy as _get_energy
from fragresp.constants import hartree_to_kcal

def check_opt(database, frag_dir, surr_cap_dir, mol_dir, stdout, stderr):

    frag_count_path = frag_dir+"/conf_count.dat"
    mol_count_path  = mol_dir+"/conf_count.dat"
    surr_count_path = surr_cap_dir+"/conf_count.dat"

    frag_check_path = frag_dir+"/check_opt.log"
    mol_check_path  = mol_dir+"/check_opt.log"
    surr_check_path = surr_cap_dir+"/check_opt.log"

    db = pickle.load(open(database, "rb"))

    include_list_mol = list()
    mol2frag         = db.get_mol2frag()
    for mol_i in range(db.get_mol_count()):
        frag_count = len(mol2frag[mol_i])
        if frag_count==0:
            include_list_mol.append(mol_i)
    
    include_list_conn   = list()
    for conn_i, conn in enumerate(db.get_conn_list()):
        if conn.get_terminal():
            continue
        include_list_conn.append(conn_i)

    frag_count_dict = OrderedDict()
    with open(frag_count_path, 'r') as f:
        for line in f:
            l = line.rstrip().lstrip().split()
            if len(l)==0:
                continue
            i,c = int(l[0]), int(l[1])
            frag_count_dict[i] = c

    _check_opt(frag_count_dict,
        mol_list=db.get_frag_list(),
        qm_dir=frag_dir,
        logfile=frag_check_path,
        stdout=stdout,
        stderr=stderr)

    plot_cycle(frag_count_dict, frag_dir)

    if len(include_list_conn)>0:
        surr_count_dict = OrderedDict()
        with open(surr_count_path, 'r') as f:
            for line in f:
                l = line.rstrip().lstrip().split()
                if len(l)==0:
                    continue
                i,c = int(l[0]), int(l[1])
                surr_count_dict[i] = c

        _check_opt(surr_count_dict,
            mol_list=db.get_surr_cap_list(),
            qm_dir=surr_cap_dir,
            logfile=surr_check_path,
            stdout=stdout,
            stderr=stderr)

        plot_cycle(surr_count_dict, surr_cap_dir)

    if len(include_list_mol)>0:
        mol_count_dict = OrderedDict()
        with open(mol_count_path, 'r') as f:
            for line in f:
                l = line.rstrip().lstrip().split()
                if len(l)==0:
                    continue
                i,c = int(l[0]), int(l[1])
                mol_count_dict[i] = c

        _check_opt(mol_count_dict,
            mol_list=db.get_mol_list(),
            qm_dir=mol_dir,
            logfile=mol_check_path,
            stdout=stdout,
            stderr=stderr)

        plot_cycle(mol_count_dict, mol_dir)


def plot_cycle(conf_dict, path):

    for i,c in conf_dict.items():
        fig_path = path\
            +"/frag%d"%i\
            +"/opt_conv.png"
        hasplot  = False
        for j in range(c) :
            ene_list = list()
            log_path     = path\
                +"/frag%d"%i\
                +"/conf%d"%j\
                +"/frag%d-conf%d_opt.log" %(i,j)
            if os.path.exists(log_path):
                _get_energy(log_path, ene_list)
                ene_list  = np.array(ene_list)
                ene_list *= hartree_to_kcal
                plt.plot(ene_list, label="Conf %d" %j)
                hasplot   = True

        if hasplot:
            plt.xlabel("Cycle")
            plt.ylabel("Energy [kcal]")
            plt.legend()
            plt.savefig(fig_path, dpi=1000)
            plt.clf()

def check_esp(database, frag_dir, surr_cap_dir, mol_dir):

    frag_count_path = frag_dir+"/conf_count.dat"
    mol_count_path  = mol_dir+"/conf_count.dat"
    surr_count_path = surr_cap_dir+"/conf_count.dat"

    frag_check_path = frag_dir+"/check_esp.log"
    mol_check_path  = mol_dir+"/check_esp.log"
    surr_check_path = surr_cap_dir+"/check_esp.log"

    db = pickle.load(open(database, "rb"))

    include_list_mol = list()
    mol2frag         = db.get_mol2frag()
    for mol_i in range(db.get_mol_count()):
        frag_count = len(mol2frag[mol_i])
        if frag_count==0:
            include_list_mol.append(mol_i)
    
    include_list_conn   = list()
    for conn_i, conn in enumerate(db.get_conn_list()):
        if conn.get_terminal():
            continue
        include_list_conn.append(conn_i)

    frag_count_dict = OrderedDict()
    with open(frag_count_path, 'r') as f:
        for line in f:
            l = line.rstrip().lstrip().split()
            if len(l)==0:
                continue
            i,c = int(l[0]), int(l[1])
            frag_count_dict[i] = c

    _check_esp(frag_count_dict,
        qm_dir=frag_dir,
        logfile=frag_check_path)

    if len(include_list_conn)>0:
        surr_count_dict = OrderedDict()
        with open(surr_count_path, 'r') as f:
            for line in f:
                l = line.rstrip().lstrip().split()
                if len(l)==0:
                    continue
                i,c = int(l[0]), int(l[1])
                surr_count_dict[i] = c

        _check_esp(surr_count_dict,
            qm_dir=surr_cap_dir,
            logfile=surr_check_path)

    if len(include_list_mol)>0:
        mol_count_dict = OrderedDict()
        with open(mol_count_path, 'r') as f:
            for line in f:
                l = line.rstrip().lstrip().split()
                if len(l)==0:
                    continue
                i,c = int(l[0]), int(l[1])
                mol_count_dict[i] = c

        _check_esp(mol_count_dict,
            qm_dir=mol_dir,
            logfile=mol_check_path)