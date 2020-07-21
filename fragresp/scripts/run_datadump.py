import pickle
import os
import time

from fragresp.utils import logger

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdmolops
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions
DrawingOptions.bondLineWidth=1.5
DrawingOptions.includeAtomNumbers=True

def datadump(database, dumpdir):

    db = pickle.load(open(database, "rb"))

    if os.path.exists(dumpdir):
        raise Warning("Caution, %s already exists. Already existing data may be overwritten.")
    else:
        os.mkdir(dumpdir)
        os.mkdir(dumpdir+"/png")

    frag2mol = db.get_frag2mol()
    frag2lcapconn = db.get_frag2lcapconn()
    frag2rcapconn = db.get_frag2rcapconn()
    mol2frag = db.get_mol2frag()
    mol2conn = db.get_mol2conn()

    frag_log = logger(dumpdir+"/frag.dat")
    frag_log.log("### datadump of database %s" %database)
    frag_log.log("### timestamp %s" %time.asctime( time.localtime(time.time()) ))
    frag_log.log("### written by run_fragresp.py datadump routine.")
    frag_log.log("###")
    frag_log.log("### ----------------- ###")
    frag_log.log("### FRAGMENT DATA LOG ###")
    frag_log.log("### ----------------- ###")
    frag_log.log("###")
    frag_log.log("# id smiles mol_id lcap_id rcap_id Natoms Nbonds Nnonhatoms Chg Nhbd Nhba Nrotbonds Nrings")

    for frag_i in range(db.get_frag_count()):
        frag    = db.get_frag(frag_i)
        Chem.SanitizeMol(frag)

        log_str = list()

        ### id
        log_str.append(str(frag_i)+" ")
        ### smiles
        log_str.append(str(Chem.MolToSmiles(frag, isomericSmiles=True))+" ")

        ### mol_id
        mol_count = len(frag2mol[frag_i])
        if mol_count==0:
            log_str.append("-1 ")
        else:
            for i in range(mol_count):
                mol_i = frag2mol[frag_i][i]
                if i<mol_count-1:
                    log_str.append(str(mol_i)+",")
                else:
                    log_str.append(str(mol_i)+" ")

        ### lcap_id
        lcap_count = len(frag2lcapconn[frag_i])
        if lcap_count==0:
            log_str.append("-1 ")
        else:
            for i in range(lcap_count):
                cap_i = frag2lcapconn[frag_i][i]
                if i<lcap_count-1:
                    log_str.append(str(cap_i)+",")
                else:
                    log_str.append(str(cap_i)+" ")

        ### rcap_id
        rcap_count = len(frag2rcapconn[frag_i])
        if rcap_count==0:
            log_str.append("-1 ")
        else:
            for i in range(rcap_count):
                cap_i = frag2rcapconn[frag_i][i]
                if i<rcap_count-1:
                    log_str.append(str(cap_i)+",")
                else:
                    log_str.append(str(cap_i)+" ")

        ### N_atoms
        log_str.append(str(frag.GetNumAtoms())+" ")
        ### N_bonds
        log_str.append(str(frag.GetNumBonds())+" ")
        ### Nnonhatoms
        log_str.append(str(frag.GetNumHeavyAtoms())+" ")
        ### Chg
        log_str.append(str(rdmolops.GetFormalCharge(frag))+" ")
        ### Nhbd
        log_str.append(str(rdMolDescriptors.CalcNumHBD(frag))+" ")
        ### Nhba
        log_str.append(str(rdMolDescriptors.CalcNumHBA(frag))+" ")
        ### Nrotbonds
        log_str.append(str(rdMolDescriptors.CalcNumRotatableBonds(frag))+" ")
        ### Nrings
        log_str.append(str(rdMolDescriptors.CalcNumRings(frag))+" ")

        frag_log.log("".join(log_str))

        png_path=dumpdir+"/png/"+"frag_%d.png" %frag_i
        try:
            Chem.SanitizeMol(frag)
            AllChem.Compute2DCoords(frag)
            Draw.MolToFile(frag,png_path, size=(500, 500))
        except:
            #Chem.Kekulize(frag)
            print ("Could not save frag %d to disk." %frag_i)

    frag_log.close()

    mol_log  = logger(dumpdir+"/mol.dat")
    mol_log.log("### datadump of database %s" %database)
    mol_log.log("### timestamp %s" %time.asctime( time.localtime(time.time()) ))
    mol_log.log("### written by run_fragresp.py datadump routine.")
    mol_log.log("###")
    mol_log.log("### ----------------- ###")
    mol_log.log("### MOLECULE DATA LOG ###")
    mol_log.log("### ----------------- ###")
    mol_log.log("###")
    mol_log.log("# id name smiles frag_id Natoms Nbonds Nnonhatoms Chg Nhbd Nhba Nrotbonds Nrings")

    for mol_i in range(db.get_mol_count()):
        mol    = db.get_mol(mol_i)
        Chem.SanitizeMol(mol)
        name   = db.get_name(mol_i)
        decomp = db.get_decompose(mol_i)

        log_str = list()

        log_str.append(str(mol_i)+" ")
        log_str.append(name+" ")
        log_str.append(str(Chem.MolToSmiles(mol, isomericSmiles=True))+" ")

        frag_count = decomp.get_frag_count()

        if frag_count==0:
            log_str.append("-1 ")
        else:
            for i in range(frag_count):
                frag_i = mol2frag[mol_i][i]
                if i<frag_count-1:
                    log_str.append(str(frag_i)+",")
                else:
                    log_str.append(str(frag_i)+" ")

        log_str.append(str(mol.GetNumAtoms())+" ")
        log_str.append(str(mol.GetNumBonds())+" ")
        log_str.append(str(mol.GetNumHeavyAtoms())+" ")
        log_str.append(str(rdmolops.GetFormalCharge(mol))+" ")
        log_str.append(str(rdMolDescriptors.CalcNumHBD(mol))+" ")
        log_str.append(str(rdMolDescriptors.CalcNumHBA(mol))+" ")
        log_str.append(str(rdMolDescriptors.CalcNumRotatableBonds(mol))+" ")
        log_str.append(str(rdMolDescriptors.CalcNumRings(mol))+" ")

        mol_log.log("".join(log_str))

        png_path=dumpdir+"/png/"+"mol_%d.png" %mol_i
        AllChem.Compute2DCoords(mol)
        Chem.Kekulize(mol)
        Draw.MolToFile(mol,png_path, size=(500, 500))

    mol_log.close()

    surr_log  = logger(dumpdir+"/surr.dat")
    surr_log.log("### datadump of database %s" %database)
    surr_log.log("### timestamp %s" %time.asctime( time.localtime(time.time()) ))
    surr_log.log("### written by run_fragresp.py datadump routine.")
    surr_log.log("###")
    surr_log.log("### ----------------- ###")
    surr_log.log("### SURROGATE DATA LOG ###")
    surr_log.log("### ------------------ ###")
    surr_log.log("###")
    surr_log.log("# id name smiles mol_id Natoms Nbonds Nnonhatoms Chg Nhbd Nhba Nrotbonds Nrings")

    for conn_i, conn in enumerate(db.get_conn_list()):

        if conn.get_terminal():
            continue

        name = conn.get_name()

        conn_cap = conn.get_surrogate_cap()
        Chem.SanitizeMol(conn_cap)

        log_str = list()

        log_str.append(str(conn_i)+" ")
        log_str.append(name+" ")
        log_str.append(str(Chem.MolToSmiles(conn_cap, isomericSmiles=True))+" ")

        conn2mol  = db.get_conn2mol()[conn_i]
        mol_count = len(conn2mol)

        if mol_count==0:
            log_str.append("-1 ")
        else:
            for i in range(mol_count):
                mol_i = conn2mol[i]
                if i<mol_count-1:
                    log_str.append(str(mol_i)+",")
                else:
                    log_str.append(str(mol_i)+" ")

        log_str.append(str(conn_cap.GetNumAtoms())+" ")
        log_str.append(str(conn_cap.GetNumBonds())+" ")
        log_str.append(str(conn_cap.GetNumHeavyAtoms())+" ")
        log_str.append(str(rdmolops.GetFormalCharge(conn_cap))+" ")
        log_str.append(str(rdMolDescriptors.CalcNumHBD(conn_cap))+" ")
        log_str.append(str(rdMolDescriptors.CalcNumHBA(conn_cap))+" ")
        log_str.append(str(rdMolDescriptors.CalcNumRotatableBonds(conn_cap))+" ")
        log_str.append(str(rdMolDescriptors.CalcNumRings(conn_cap))+" ")

        surr_log.log("".join(log_str))

        png_path=dumpdir+"/png/"+"surr_%s.png" %(conn_i)
        AllChem.Compute2DCoords(conn_cap)
        Chem.Kekulize(conn_cap)
        Draw.MolToFile(conn_cap, png_path, size=(500, 500))

    surr_log.close()