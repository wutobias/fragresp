import pickle
from collections import OrderedDict

from fragresp.resp_fit import resp_it_db

def resp(database, 
    frag_dir,
    surr_cap_dir,
    mol_dir,
    fit_to_surr,
    remap,
    by_charge,
    capsum,
    remap_dir,
    stdout,
    stderr):

    db = pickle.load(open(database, "rb"))

    frag_count_path = frag_dir+"/conf_count.dat"
    mol_count_path  = mol_dir+"/conf_count.dat"
    surr_count_path = surr_cap_dir+"/conf_count.dat"

    frag_count_dict   = OrderedDict()
    mol_count_dict    = OrderedDict()
    surr_count_dict   = OrderedDict()

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

    with open(frag_count_path, 'r') as f:
        for line in f:
            l = line.rstrip().lstrip().split()
            if len(l)==0:
                continue
            i,c = int(l[0]), int(l[1])
            frag_count_dict[i] = c

    if len(include_list_conn)>0:
        with open(surr_count_path, 'r') as f:
            for line in f:
                l = line.rstrip().lstrip().split()
                if len(l)==0:
                    continue
                i,c = int(l[0]), int(l[1])
                surr_count_dict[i] = c

    if len(include_list_mol)>0:
        with open(mol_count_path, 'r') as f:
            for line in f:
                l = line.rstrip().lstrip().split()
                if len(l)==0:
                    continue
                i,c = int(l[0]), int(l[1])
                mol_count_dict[i] = c

    resp_it_db(db, frag_dir, frag_count_dict, surr_cap_dir, surr_count_dict,
        mol_dir, mol_count_dict, fit_to_surr, remap, by_charge, capsum, 
        remap_dir, stdout, stderr)