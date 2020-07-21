import pickle
from collections import OrderedDict

def conformers(database, 
                add_db, 
                limit, 
                percentage, 
                frag_dir, 
                surr_cap_dir, 
                mol_dir, 
                fraglog, 
                surrlog, 
                mollog, 
                pipeline,
                queue):

    if pipeline=='openbabel':
        from fragresp.babel_conformers import prep_qm
    else:
        from fragresp.oe_conformers import prep_qm

    overwrite = not add_db

    frag_count_path = frag_dir+"/conf_count.dat"
    surr_count_path = surr_cap_dir+"/conf_count.dat"
    mol_count_path  = mol_dir+"/conf_count.dat"

    db = pickle.load(open(database, "rb"))

    include_list_mol = list()
    mol2frag         = db.get_mol2frag()
    for mol_i in range(db.get_mol_count()):
        frag_count = len(mol2frag[mol_i])
        if frag_count==0:
            include_list_mol.append(mol_i)
    
    include_list_conn   = list()
    for conn_i, conn in enumerate(db.get_conn_list()):
        if not conn.get_terminal():
            include_list_conn.append(conn_i)

    frag_count = prep_qm(db.get_frag_list(), 
                        overwrite, 
                        limit, 
                        percentage,
                        qm_dir=frag_dir, 
                        logfile=fraglog,
                        queue=queue)

    if len(include_list_conn)>0:
        surr_count = prep_qm(db.get_surr_cap_list(),
                            overwrite, 
                            limit, 
                            percentage,
                            qm_dir=surr_cap_dir, 
                            logfile=surrlog,
                            include_list=include_list_conn,
                            queue=queue)
    
    if len(include_list_mol)>0:
        mol_count = prep_qm(db.get_mol_list(),
                            overwrite, 
                            limit, 
                            percentage,
                            qm_dir=mol_dir, 
                            logfile=mollog,
                            include_list=include_list_mol,
                            queue=queue)

    if overwrite:
        with open(frag_count_path, 'w') as f:
            for i, c in enumerate(frag_count):
                f.write("%d %d\n" %(i,c))

        if len(include_list_conn)>0:
            with open(surr_count_path, 'w') as f:
                for i, c in zip(include_list_conn, surr_count):
                    f.write("%d %d\n" %(i,c))

        if len(include_list_mol)>0:
            with open(mol_count_path, 'w') as f:
                for i, c in zip(include_list_mol, mol_count):
                    f.write("%d %d\n" %(i,c))

    else:
        ### FragCount file ###
        ### -------------- ###
        frag_count_old = OrderedDict()
        with open(frag_count_path, 'r') as f:
            for line in f:
                l = line.rstrip().lstrip().split()
                if len(l)==0:
                    continue
                i,c = int(l[0]), int(l[1])
                frag_count_old[i] = c
        new_frag_str = list()
        for i, c in enumerate(frag_count):
            ### if c equals -1, then no new conformer
            ### was generated for this fragment. In this
            ### case we take the value from the old frag_count file.
            if c == -1:
                if i not in frag_count_old:
                    raise Warning("Conformer count for fragment %d not found in %s." %(i,frag_count_path))
                new_frag_str.append("%d %d\n" %(i,frag_count_old[i]))
            else:
                new_frag_str.append("%d %d\n" %(i,c))

        with open(frag_count_path, 'w') as f:
            f.write("".join(new_frag_str))
        
        if len(include_list_conn)>0:
            ### SurrCount file ###
            ### -------------- ###
            surr_count_old = OrderedDict()
            with open(surr_count_path, 'r') as f:
                for line in f:
                    l = line.rstrip().lstrip().split()
                    if len(l)==0:
                        continue
                    i,c = int(l[0]), int(l[1])
                    surr_count_old[i] = c
            new_surr_str = list()
            for i, c in zip(include_list_conn, surr_count):
                ### if c equals -1, then no new conformer
                ### was generated for this fragment. In this
                ### case we take the value from the old surr_count file.
                if c == -1:
                    if i not in surr_count_old:
                        raise Warning("Conformer count for surrogate %d not found in %s." %(i,surr_count_path))
                    new_surr_str.append("%d %d\n" %(i,surr_count_old[i]))
                else:
                    new_surr_str.append("%d %d\n" %(i,c))

            with open(surr_count_path, 'w') as f:
                f.write("".join(new_surr_str))

        if len(include_list_mol)>0:
            ### MolCount file ###
            ### ------------- ###
            mol_count_old = OrderedDict()
            with open(mol_count_path, 'r') as f:
                for line in f:
                    l = line.rstrip().lstrip().split()
                    if len(l)==0:
                        continue
                    i,c = int(l[0]), int(l[1])
                    mol_count_old[i] = c
            new_mol_str = list()
            for i, c in zip(include_list_mol, mol_count):
                ### if c equals -1, then no new conformer
                ### was generated for this fragment. In this
                ### case we take the value from the old surr_count file.
                if c == -1:
                    if i not in surr_count_old:
                        raise Warning("Conformer count for molecule %d not found in %s." %(i,mol_count_path))
                    new_mol_str.append("%d %d\n" %(i,surr_count_old[i]))
                else:
                    new_mol_str.append("%d %d\n" %(i,c))

            with open(mol_count_path, 'w') as f:
                f.write("".join(new_mol_str))