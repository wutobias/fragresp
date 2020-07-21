import pickle

from fragresp.decomposition import decompose
from fragresp.data_strc import fragment_db

from rdkit import Chem

def decompose_molecules(mol2_list, database_path, linker_list, add_db, verbose):

    if not add_db:
        db = fragment_db(verbose)
    else:
        ### Load already existing db
        db = pickle.load(open(database_path, "rb"))

    for i, mol2 in enumerate(mol2_list):
        if verbose:
            print ("mol2list %d:" %i)
        j =0
        for name, path in mol2.items():

            path = path[0]

            if verbose:
                print ("decompose %s %s" %(name, path))

            mol          = Chem.MolFromMol2File(path, removeHs=True)
            decompose_it = decompose(mol, verbose)

            if add_db:
                for linker in db.get_conn_list():
                    decompose_it.add_connector(linker)

            else:
                for linker in linker_list:
                    decompose_it.add_connector(linker)

            decompose_it.decompose()

            frag_count = db.get_frag_count()
            if verbose:
                print ("Attempting addition to database...")
            if db.add_mol(decompose_it, path, name):
                if add_db or verbose:
                    print ("--database:")
                    print ("Molecule %s new to database." %name)
                    if frag_count != db.get_frag_count():
                        print ("Fragments",)
                        for frag_id in range(frag_count, db.get_frag_count()):
                            print (frag_id,)
                        print ("new to database.")

    pickle.dump(db, open(database_path, "wb"))