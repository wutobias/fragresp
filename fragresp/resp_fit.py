import os
from subprocess import call
import fortranformat
import numpy as np
from rdkit import Chem

from fragresp.aux_progs import resp_exe
from fragresp.utils import mol2_reader
from fragresp.utils import cat_esp
from fragresp.utils import preprocess_esp
from fragresp.utils import fix_groups
from fragresp.resp_utils import write_resp_in
from fragresp.resp_utils import resp_surrogate
from fragresp.gaussian_tools import get_energy
from fragresp.constants import hartree_to_kcal

### Suppress warnings when loading mol2 files without
### hydrogen atoms into rdkit.
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

def resp_it_db(db,
            frag_path,
            frag_conf_dict, 
            surr_path=None,
            surr_conf_dict=None,
            mol_path=None,
            mol_conf_dict=None,
            fit_to_surr=False,
            remap=True,
            by_charge=True,
            capsum=True,
            remap_dir="remap_dir",
            energy_weighting=False,
            stdout=None,
            stderr=None):

    if fit_to_surr:
        if type(surr_path) == None:
            raise Exception("If fit_to_surr=True, surr_path must be given.")
        if type(surr_conf_dict) == None:
            raise Exception("If fit_to_surr=True, surr_conf_dict must be given.")

    if stdout==None:
        stdout=open(os.devnull, 'w')
    if stderr==None:
        stderr=open(os.devnull, 'w')

    ####################################################
    ###                                              ###
    ### First, calculate resp charges for the capped ###
    ### surrogate structure.                         ###
    ###                                              ###
    ####################################################

    N_conn = db.get_conn_count()

    surr_mol2path_dict    = dict()
    surr_mol_dict         = dict()
    surr_mol_matches_dict = dict()

    if fit_to_surr:

        for conn_i in range(N_conn):

            conn_surr = db.get_conn(conn_i)
            if conn_surr.get_terminal():
                continue
            surr_path_main = surr_path+"/frag%d" %conn_i
            mol2_out_path, mol_qm_surr, matches = resp_surrogate(conn_surr,
                                                                surr_conf_dict,
                                                                surr_path_main,
                                                                conn_i,
                                                                energy_weighting,
                                                                stdout,
                                                                stderr,
                                                                "None",
                                                                capsum)
            surr_mol2path_dict[conn_i]    = mol2_out_path
            surr_mol_dict[conn_i]         = mol_qm_surr
            surr_mol_matches_dict[conn_i] = matches


    ########################################################
    ###                                                  ###
    ### Second, calculate resp charges for the capped    ###
    ### fragments, that are found for the each connector ###
    ###                                                  ###
    ########################################################

    N_frag = db.get_frag_count()

    ### chg_list   : store all different charges found in the
    ###              fragment database.
    chg_list = list()

    if by_charge:
        for frag_i in range(N_frag):

            ### Cluster the whole fragment database into
            ### sets of equal charge.
            mol_frag   = db.get_frag(frag_i)
            charge     = Chem.GetFormalCharge(mol_frag)
            if charge not in chg_list:
                chg_list.append(charge)

    else:
        chg_list.append(0)

    for chg_i in range(len(chg_list)):

        if by_charge:
            charge = chg_list[chg_i]

        ### conn_id_rcap_list      : store connector idx of some rcap
        ### conn_rcap_id_list      : store fragment id of first occurence
        ###                          of a particular rcap in some connector id
        ### conn_rcap_atm_idxs_list: store atom indices of rcap in this fragment
        ### conn_ranc_idxs_list      : store atom indices of anchor atoms that have frozen
        ###                          charges in that fragment.

        conn_id_rcap_list  = list()
        conn_rcap_id_list  = list()
        conn_rcap_atm_idxs_list = list()
        conn_ranc_idxs_list  = list()
        conn_ranc_idxs_list_noh = list()

        ### conn_id_lcap_list      : store connector idx of some lcap
        ### conn_lcap_id_list      : store fragment id of first occurence
        ###                          of a particular lcap in some connector id
        ### conn_lcap_atm_idxs_list: store atom indices of rcap in this fragment
        ### conn_lanc_idxs_list      : store atom indices of anchor atoms that have frozen
        ###                            charges in that fragment.
        ### fragdb2fragresp       : Stores the entry of the fragment in the resp input (last
        ###                         conformer) at the position of the entry in the frag_db
        ### fragdb2fragresp_matches : Atom matches for atoms in fragdb and resp input
        ### fragdb2fragresp_mol     : mol objects as obtained from QM calculation for each
        ###                           entry in the fragdb

        conn_id_lcap_list  = list()
        conn_lcap_id_list  = list()
        conn_lcap_atm_idxs_list = list()
        conn_lanc_idxs_list  = list()
        conn_lanc_idxs_list_noh = list()

        fragdb2fragresp = list()
        fragdb2fragresp_matches = list()
        fragdb2fragresp_mol     = list()

        ### The input files for the resp program will be generated
        ### with the write_resp_in file generator.
        frag_resp_in = write_resp_in()
        esp_list     = list()

        if by_charge:
            charge_path_main  = frag_path+"/chg%s" %charge
        else:
            charge_path_main  = frag_path+"/chgX"

        if not os.path.exists(charge_path_main):
            os.mkdir(charge_path_main)

        charge_path_esp   = charge_path_main+"/esp"
        charge_path_resp1 = charge_path_main+"/resp1"
        charge_path_resp2 = charge_path_main+"/resp2"

        frag_i_int = 0

        for frag_i in range(N_frag):

            frag_path_main  = frag_path+"/frag%d" %frag_i
            frag_path_esp   = frag_path_main+"/esp"

            mol_frag   = db.get_frag(frag_i)

            if by_charge:
                ### Only fragments with equal charges
                if charge != Chem.GetFormalCharge(mol_frag):
                    continue

            conformer_weights = list()
            if energy_weighting:
                for conf_i in range(frag_conf_dict[frag_i]):
                    conf_i_path = frag_path_main+"/conf%d/frag%d-conf%d_esp.log" %(conf_i, frag_i, conf_i)
                    energy_list         = list()
                    get_energy(conf_i_path, energy_list)
                    conformer_weights.append(
                        energy_list[-1] * hartree_to_kcal
                        )
                min_energy = min(conformer_weights)
                for idx in range(len(conformer_weights)):
                    conformer_weights[idx] -= min_energy
                    ### gas constant in cal/mol/K
                    conformer_weights[idx] = np.exp(-conformer_weights[idx]/(1.9872159 * 1.e-3 * 300.))
            else:
                for conf_i in range(frag_conf_dict[frag_i]):
                    conformer_weights.append(1.)

            frag_i_ref = frag_i_int
            for idx, conf_i in enumerate(range(frag_conf_dict[frag_i])):
                conf_i_path = frag_path_main+"/conf%d/frag%d-conf%d_esp.log" %(conf_i, frag_i, conf_i)
                esp_path, mol2_path = preprocess_esp(conf_i_path)
                ### Match the mol2 structure obtained from the QM calculation
                ### with the rdkit mol object stored in the database. Perform
                ### some sanity checks on the match lists.
                mol_qm_frag = Chem.MolFromMol2File(mol2_path, removeHs=False)
                mol_qm_frag = fix_groups(mol_qm_frag)
                mol_frag    = fix_groups(mol_frag)
                matches     = mol_qm_frag.GetSubstructMatches(mol_frag)

                if len(matches) == 0:
                    raise Exception("Did not find any substructure matches.")
                if len(matches) > 1:
                    raise Exception("Found %d substructure matches. Should be one only." %(len(matches)))
                if len(matches[0]) != mol_frag.GetNumAtoms():
                    raise Exception("Found %d atoms in substructure match, but expected %d." %(len(matches[0]),\
                                     mol_frag.GetNumAtoms()))
                esp_list.append(esp_path)
                frag_resp_in.add_mol(mol_qm_frag,\
                                     Chem.GetFormalCharge(mol_qm_frag),\
                                     conformer_weights[idx],\
                                     "Frag %d Conf %d, resp mol %d" %(frag_i,\
                                                                      conf_i,\
                                                                      frag_resp_in._mol_count+1))
                ### Generate intermolecular restraints
                frag_resp_in.add_intermolecular(frag_i_ref, frag_i_int)

                frag_i_int += 1

            matches = matches[0]

            fragdb2fragresp.append(frag_i_ref)
            fragdb2fragresp_matches.append(matches)
            fragdb2fragresp_mol.append(mol_qm_frag)

            ### Set atom center charge restraints (freeze atoms)
            ### We will freeze only those atoms, that are part of the
            ### anchor (ranc or lanc) of a fragment and are not 
            ### explicitly set vary freely.

            for cap_conn, cap_idxs, anc_idxs in zip(db.get_frag2rcapconn()[frag_i],\
                                                    db.get_frag2rcapmap()[frag_i],\
                                                    db.get_frag2rancmap()[frag_i]):

                conn = db.get_conn(cap_conn)

                conn_anc_map = conn.get_ranc_map()
                if conn.get_terminal():
                    conn_cap_map = list()
                else:
                    conn_cap_map = conn.get_rcap_map()

                #print ("R site")
                #print ("frag_i", frag_i)
                #print (conn.get_name())
                #if conn.get_terminal():
                #    print ("I'm terminal dude.")
                #print ("cap_conn", cap_conn)
                #print ("cap_idxs", cap_idxs)
                #print ("anc_idxs", anc_idxs)

                ### R anchor and R capping groups
                ### -----------------------------
                if cap_conn not in conn_id_rcap_list:

                    conn_id_rcap_list.append(cap_conn)
                    conn_rcap_id_list.append(frag_i_ref)
                    conn_ranc_idxs_list.append(list())
                    conn_ranc_idxs_list_noh.append(list())
                    conn_rcap_atm_idxs_list.append(list())

                    if fit_to_surr and not conn.get_terminal():
                        surr_mol2    = mol2_reader(surr_mol2path_dict[cap_conn])
                        surr_mol     = surr_mol_dict[cap_conn]
                        surr_matches = surr_mol_matches_dict[cap_conn]

                    ### Make charges on anchor atoms equal in all
                    ### molecules that are bound to the same cap
                    ### at the same connector.
                    for anc_idx, anc_conn in zip(anc_idxs,
                                                conn_anc_map):
                        if anc_conn < 0:
                            continue
                        frag_atm_idx = matches[anc_idx]
                        atm          = mol_qm_frag.GetAtomWithIdx(frag_atm_idx)
                        conn_ranc_idxs_list[-1].append(frag_atm_idx)
                        conn_ranc_idxs_list_noh[-1].append(frag_atm_idx)
                        for neighbor in atm.GetNeighbors():
                            neighbor_idx = neighbor.GetIdx()
                            if neighbor.GetAtomicNum() == 1 \
                            and neighbor_idx not in conn_ranc_idxs_list[-1]:
                                conn_ranc_idxs_list[-1].append(neighbor_idx)

                    ### Add restraints for whole capping group
                    atom_i_list   = list()
                    mol_i_list    = list()
                    cap_neighbors = list()
                    for cap_idx in cap_idxs:
                        cap_neighbors.append(matches[cap_idx])

                    if fit_to_surr and not conn.get_terminal():
                        cap_neighbors_surr = list()
                        for cap_idx in conn_cap_map:
                            cap_neighbors_surr.append(surr_matches[cap_idx])

                    for cap_idx, conn_cap_idx in zip(cap_idxs,
                                                     conn_cap_map):
                        atm_idx = matches[cap_idx]
                        atm     = mol_qm_frag.GetAtomWithIdx(atm_idx)

                        if fit_to_surr and not conn.get_terminal():
                            surr_cap_idx = surr_matches[conn_cap_idx]
                            atm_surr     = surr_mol.GetAtomWithIdx(surr_cap_idx)
                            frag_resp_in.add_group([atm_idx],
                                                   [frag_i_ref],
                                                    surr_mol2.charges[surr_cap_idx])
                            frag_resp_in.add_free(frag_i_ref, [atm_idx])

                        conn_rcap_atm_idxs_list[-1].append(atm_idx)
                        mol_i_list.append(frag_i_ref)
                        matched_neighbors = list()
                        for neighbor in atm.GetNeighbors():
                            neighbor_idx = neighbor.GetIdx()
                            if neighbor.GetAtomicNum() == 1 \
                            and neighbor_idx not in cap_neighbors:
                                
                                if not fit_to_surr and not conn.get_terminal():
                                    conn_rcap_atm_idxs_list[-1].append(neighbor_idx)
                                    mol_i_list.append(frag_i_ref)

                                else:
                                    ### Here we search for corresponding hydrogen atom in the
                                    ### surrogate molecule. We do not care about which hydrogen
                                    ### atom we choose from the surrogate molecule. However,
                                    ### we take care that each hydrogen atom is choosen only once.
                                    for neighbor_surr in atm_surr.GetNeighbors():
                                        neighbor_surr_idx = neighbor_surr.GetIdx()
                                        if neighbor_surr.GetAtomicNum() == 1 \
                                        and neighbor_surr_idx not in cap_neighbors_surr \
                                        and neighbor_surr_idx not in matched_neighbors:
                                            conn_rcap_atm_idxs_list[-1].append(neighbor_idx)
                                            mol_i_list.append(frag_i_ref)
                                            frag_resp_in.add_group([neighbor_idx],
                                                                   [frag_i_ref],
                                                                   surr_mol2.charges[neighbor_surr_idx])
                                            frag_resp_in.add_free(frag_i_ref, [neighbor_idx])
                                            matched_neighbors.append(neighbor_surr_idx)
                                            break

                    if not fit_to_surr \
                    and not capsum \
                    and not conn.get_terminal():
                        if len(conn_rcap_atm_idxs_list[-1])>1:
                            frag_resp_in.add_group(conn_rcap_atm_idxs_list[-1],
                                mol_i_list,
                                0.0)
                    for idx, mol_i in enumerate(mol_i_list):
                        frag_resp_in.add_free(mol_i, [conn_rcap_atm_idxs_list[-1][idx]])

                else:

                    ### All other fragments are referenced to the 
                    ### definition of the initial restraints for this atomic center
                    index          = conn_id_rcap_list.index(cap_conn)
                    rcap_conn_id   = conn_rcap_id_list[index]
                    rcap_atm_idxs  = conn_rcap_atm_idxs_list[index]
                    ranc_atom_idxs = conn_ranc_idxs_list[index]
                    ranc_atom_idxs_noh = conn_ranc_idxs_list_noh[index]

                    atm_idx_list  = list()
                    atm_idx_list_noh = list()
                    ### Anchor atoms...
                    for anc_idx, anc_conn in zip(anc_idxs,
                                                 conn_anc_map):
                        if anc_conn < 0:
                            continue

                        frag_atm_idx = matches[anc_idx]
                        atm          = mol_qm_frag.GetAtomWithIdx(frag_atm_idx)
                        atm_idx_list.append(frag_atm_idx)
                        atm_idx_list_noh.append(frag_atm_idx)
                        for neighbor in atm.GetNeighbors():
                            neighbor_idx = neighbor.GetIdx()
                            if neighbor.GetAtomicNum() == 1 \
                            and neighbor_idx not in atm_idx_list:
                                atm_idx_list.append(neighbor_idx)

                    if len(ranc_atom_idxs) != len(atm_idx_list):
                        frag_resp_in.add_atom_equiv(ranc_atom_idxs_noh,
                                                    atm_idx_list_noh,
                                                    [rcap_conn_id]*len(ranc_atom_idxs_noh),
                                                    [frag_i_ref]*len(atm_idx_list_noh))

                    elif len(ranc_atom_idxs)>0:
                        frag_resp_in.add_atom_equiv(ranc_atom_idxs,
                                                    atm_idx_list,
                                                    [rcap_conn_id]*len(ranc_atom_idxs),
                                                    [frag_i_ref]*len(atm_idx_list))

                    ### Capping groups...
                    atm_idx_list  = list()
                    cap_neighbors = list()

                    for cap_idx in cap_idxs:
                        cap_neighbors.append(matches[cap_idx])
                    for cap_idx in cap_idxs:
                        atm_idx = matches[cap_idx]
                        atm     = mol_qm_frag.GetAtomWithIdx(atm_idx)
                        atm_idx_list.append(atm_idx)
                        for neighbor in atm.GetNeighbors():
                            neighbor_idx = neighbor.GetIdx()
                            if neighbor.GetAtomicNum() == 1 \
                               and not neighbor_idx in cap_neighbors:
                                atm_idx_list.append(neighbor_idx)

                    if len(rcap_atm_idxs)>0:
                        frag_resp_in.add_atom_equiv(rcap_atm_idxs,
                                                    atm_idx_list,
                                                    [rcap_conn_id]*len(rcap_atm_idxs),
                                                    [frag_i_ref]*len(atm_idx_list))

            ### End Loop over Rcap

            for cap_conn, cap_idxs, anc_idxs in zip(db.get_frag2lcapconn()[frag_i],\
                                                    db.get_frag2lcapmap()[frag_i],\
                                                    db.get_frag2lancmap()[frag_i]):

                conn = db.get_conn(cap_conn)

                conn_anc_map = conn.get_lanc_map()
                if conn.get_terminal():
                    conn_cap_map = list()
                else:
                    conn_cap_map = conn.get_lcap_map()

                #print ("L site")
                #print ("frag_i", frag_i)
                #print (conn.get_name())
                #if conn.get_terminal():
                #    print ("I'm terminal dude.")
                #print ("cap_conn", cap_conn)
                #print ("cap_idxs", cap_idxs)
                #print ("anc_idxs", anc_idxs)

                ### L anchor and L capping groups
                ### -----------------------------
                if cap_conn not in conn_id_lcap_list:

                    conn_id_lcap_list.append(cap_conn)
                    conn_lcap_id_list.append(frag_i_ref)
                    conn_lanc_idxs_list.append(list())
                    conn_lanc_idxs_list_noh.append(list())
                    conn_lcap_atm_idxs_list.append(list())

                    if fit_to_surr and not conn.get_terminal():
                        surr_mol2    = mol2_reader(surr_mol2path_dict[cap_conn])
                        surr_mol     = surr_mol_dict[cap_conn]
                        surr_matches = surr_mol_matches_dict[cap_conn]

                    ### Make charges on anchor atoms equal in all
                    ### molecules that are bound to the same cap
                    ### at the same connector.
                    for anc_idx, anc_conn in zip(anc_idxs,
                                                 conn_anc_map):
                        if anc_conn < 0:
                            continue
                        frag_atm_idx = matches[anc_idx]
                        atm          = mol_qm_frag.GetAtomWithIdx(frag_atm_idx)
                        conn_lanc_idxs_list[-1].append(frag_atm_idx)
                        conn_lanc_idxs_list_noh[-1].append(frag_atm_idx)
                        for neighbor in atm.GetNeighbors():
                            neighbor_idx = neighbor.GetIdx()
                            if neighbor.GetAtomicNum() == 1 \
                            and neighbor_idx not in conn_lanc_idxs_list[-1]:
                                conn_lanc_idxs_list[-1].append(neighbor_idx)

                    ### Add restraints for whole capping group
                    atom_i_list = list()
                    mol_i_list  = list()
                    cap_neighbors = list()
                    for cap_idx in cap_idxs:
                        cap_neighbors.append(matches[cap_idx])

                    if fit_to_surr and not conn.get_terminal():
                        cap_neighbors_surr = list()
                        for cap_idx in conn_cap_map:
                            cap_neighbors_surr.append(surr_matches[cap_idx])

                    for cap_idx, conn_cap_idx in zip(cap_idxs,
                                                     conn_cap_map):
                        atm_idx = matches[cap_idx]
                        atm     = mol_qm_frag.GetAtomWithIdx(atm_idx)
                        
                        if fit_to_surr and not conn.get_terminal():
                            surr_cap_idx = surr_matches[conn_cap_idx]
                            atm_surr     = surr_mol.GetAtomWithIdx(surr_cap_idx)
                            frag_resp_in.add_group([atm_idx],
                                                   [frag_i_ref],
                                                    surr_mol2.charges[surr_cap_idx])
                            frag_resp_in.add_free(frag_i_ref, [atm_idx])

                        conn_lcap_atm_idxs_list[-1].append(atm_idx)
                        mol_i_list.append(frag_i_ref)
                        matched_neighbors = list()
                        for neighbor in atm.GetNeighbors():
                            neighbor_idx = neighbor.GetIdx()
                            if neighbor.GetAtomicNum() == 1 \
                            and neighbor_idx not in cap_neighbors:
                                
                                if not fit_to_surr and not conn.get_terminal():
                                    conn_lcap_atm_idxs_list[-1].append(neighbor_idx)
                                    mol_i_list.append(frag_i_ref)

                                else:
                                    ### Here we search for corresponding hydrogen atom in the
                                    ### surrogate molecule. We do not care about which hydrogen
                                    ### atom we choose from the surrogate molecule. However,
                                    ### we take care that each hydrogen atom is choosen only once.
                                    for neighbor_surr in atm_surr.GetNeighbors():
                                        neighbor_surr_idx = neighbor_surr.GetIdx()
                                        if neighbor_surr.GetAtomicNum() == 1 \
                                        and neighbor_surr_idx not in cap_neighbors_surr \
                                        and neighbor_surr_idx not in matched_neighbors:
                                            conn_lcap_atm_idxs_list[-1].append(neighbor_idx)
                                            mol_i_list.append(frag_i_ref)
                                            frag_resp_in.add_group([neighbor_idx], 
                                                                    [frag_i_ref],
                                                                    surr_mol2.charges[neighbor_surr_idx])
                                            frag_resp_in.add_free(frag_i_ref, [neighbor_idx])
                                            matched_neighbors.append(neighbor_surr_idx)
                                            break

                    if not fit_to_surr \
                    and not capsum \
                    and not conn.get_terminal():
                        if len(conn_lcap_atm_idxs_list[-1])>1:
                            frag_resp_in.add_group(conn_lcap_atm_idxs_list[-1],
                                mol_i_list,
                                0.0)
                    for idx, mol_i in enumerate(mol_i_list):
                        frag_resp_in.add_free(mol_i, [conn_lcap_atm_idxs_list[-1][idx]])

                else:
                    ### All other fragments just refer to the
                    ### definition of the initial restraints for this atomic center
                    index          = conn_id_lcap_list.index(cap_conn)
                    lcap_conn_id   = conn_lcap_id_list[index]
                    lcap_atm_idxs  = conn_lcap_atm_idxs_list[index]
                    lanc_atom_idxs = conn_lanc_idxs_list[index]
                    lanc_atom_idxs_noh = conn_lanc_idxs_list_noh[index]

                    atm_idx_list  = list()
                    atm_idx_list_noh = list()
                    ### Anchor atoms...
                    for anc_idx, anc_conn in zip(anc_idxs,
                                                 conn_anc_map):
                        if anc_conn < 0:
                            continue

                        frag_atm_idx = matches[anc_idx]
                        atm          = mol_qm_frag.GetAtomWithIdx(frag_atm_idx)
                        atm_idx_list.append(frag_atm_idx)
                        for neighbor in atm.GetNeighbors():
                            neighbor_idx = neighbor.GetIdx()
                            if neighbor.GetAtomicNum() == 1 \
                            and neighbor_idx not in atm_idx_list:
                                atm_idx_list.append(neighbor_idx)

                    if len(lanc_atom_idxs) != len(atm_idx_list):
                        frag_resp_in.add_atom_equiv(lanc_atom_idxs_noh,
                                                    atm_idx_list_noh,
                                                    [lcap_conn_id]*len(lanc_atom_idxs_noh),
                                                    [frag_i_ref]*len(atm_idx_list_noh))

                    elif len(lanc_atom_idxs)>0:
                        frag_resp_in.add_atom_equiv(lanc_atom_idxs,
                                                    atm_idx_list,
                                                    [lcap_conn_id]*len(lanc_atom_idxs),
                                                    [frag_i_ref]*len(atm_idx_list))

                    ### Capping groups...
                    atm_idx_list  = list()
                    cap_neighbors = list()
                    for cap_idx in cap_idxs:
                        cap_neighbors.append(matches[cap_idx])
                    for cap_idx in cap_idxs:
                        atm_idx = matches[cap_idx]
                        atm     = mol_qm_frag.GetAtomWithIdx(atm_idx)
                        atm_idx_list.append(atm_idx)
                        for neighbor in atm.GetNeighbors():
                            neighbor_idx = neighbor.GetIdx()
                            if neighbor.GetAtomicNum() == 1 \
                               and neighbor_idx not in cap_neighbors:
                                atm_idx_list.append(neighbor_idx)

                    if len(lcap_atm_idxs)>0:
                        frag_resp_in.add_atom_equiv(lcap_atm_idxs,
                                                    atm_idx_list,
                                                    [lcap_conn_id]*len(lcap_atm_idxs),
                                                    [frag_i_ref]*len(atm_idx_list))

            ### End Loop over Lcap

        ### Add cross fragment restraints
        frag2frag_frgs1 = list()
        frag2frag_frgs2 = list()
        frag2frag_atms1 = list()
        frag2frag_atms2 = list()
        for frag2frag_frgs, frag2frag_atms in zip(db.get_frag2frag_frgs(),
                                                    db.get_frag2frag_atms()):

            mol_qm_frag1 = fragdb2fragresp_mol[frag2frag_frgs[0]]
            mol_qm_frag2 = fragdb2fragresp_mol[frag2frag_frgs[1]]
            matches1     = fragdb2fragresp_matches[frag2frag_frgs[0]]
            matches2     = fragdb2fragresp_matches[frag2frag_frgs[1]]

            frag2frag_frgs1.append(fragdb2fragresp[frag2frag_frgs[0]])
            frag2frag_frgs2.append(fragdb2fragresp[frag2frag_frgs[1]])
            frag2frag_atms1.append(matches1[frag2frag_atms[0]])
            frag2frag_atms2.append(matches2[frag2frag_atms[1]])

            atm1 = mol_qm_frag1.GetAtomWithIdx(frag2frag_atms1[-1])
            atm2 = mol_qm_frag2.GetAtomWithIdx(frag2frag_atms2[-1])

            for neighbor1, neighbor2 in zip(atm1.GetNeighbors(),
                                            atm2.GetNeighbors()):
                neighbor1_idx = neighbor1.GetIdx()
                neighbor2_idx = neighbor2.GetIdx()
                if neighbor1.GetAtomicNum() == 1 \
                and neighbor2.GetAtomicNum() == 1:
                    frag2frag_frgs1.append(fragdb2fragresp[frag2frag_frgs[0]])
                    frag2frag_frgs2.append(fragdb2fragresp[frag2frag_frgs[1]])
                    frag2frag_atms1.append(neighbor1_idx)
                    frag2frag_atms2.append(neighbor2_idx)

        frag_resp_in.add_atom_equiv(frag2frag_atms1,
                                    frag2frag_atms2,
                                    frag2frag_frgs1,
                                    frag2frag_frgs2)


        if capsum:
            for cap_i in range(len(conn_lcap_atm_idxs_list)):
                _atm_idx_list = list()
                _mol_idx_list = list()
                if conn_lcap_atm_idxs_list:
                    _atm_idx_list.extend(conn_lcap_atm_idxs_list[cap_i])
                    _mol_idx_list.extend(
                        [conn_lcap_id_list[cap_i]]*len(conn_lcap_atm_idxs_list[cap_i])
                        )                    
                if conn_rcap_atm_idxs_list:
                    _atm_idx_list.extend(conn_rcap_atm_idxs_list[cap_i])
                    _mol_idx_list.extend(
                        [conn_rcap_id_list[cap_i]]*len(conn_rcap_atm_idxs_list[cap_i])
                        )

                frag_resp_in.add_group(_atm_idx_list,
                                       _mol_idx_list,
                                       0.0)
                if conn_lcap_id_list:
                    frag_resp_in.add_free(
                        conn_lcap_id_list[cap_i], 
                        conn_lcap_atm_idxs_list[cap_i]
                        )
                if conn_rcap_id_list:
                    frag_resp_in.add_free(
                        conn_rcap_id_list[cap_i], 
                        conn_rcap_atm_idxs_list[cap_i]
                        )

            #print (frag_i_ref+1, "-->", frag_i)

        cat_esp(esp_list, charge_path_esp)
        frag_resp_in.set_title('Stage 1 of multimolecule multiconformation resp fit.')
        frag_resp_in.set_stage('resp1')
        frag_resp_in.write(charge_path_resp1)

        frag_resp_in.set_title('Stage 2 of multimolecule multiconformation resp fit.')
        frag_resp_in.set_stage('resp2')
        if fit_to_surr:
            frag_resp_in.set_stage('resp2')
        else:
            frag_resp_in.set_stage('resp2_surr')
        frag_resp_in.write(charge_path_resp2)

        resp_args = [resp_exe, "-O", "-i", charge_path_resp1, 
                               "-o", charge_path_main+"/resp1.out",
                               "-e", charge_path_esp,
                               "-t", charge_path_main+"/resp1.crg", 
                               "-p", charge_path_main+"/punch1",
                               "-s", charge_path_main+"/esout1"]
        print(" ".join(resp_args))
        call(resp_args, stdout=stdout, stderr=stderr)

        resp_args = [resp_exe, "-O", "-i", charge_path_resp2, 
                               "-o", charge_path_main+"/resp2.out",
                               "-q", charge_path_main+"/resp1.crg",
                               "-e", charge_path_esp, 
                               "-t", charge_path_main+"/resp2.crg", 
                               "-p", charge_path_main+"/punch2",
                               "-s", charge_path_main+"/esout2"]
        print(" ".join(resp_args))
        call(resp_args, stdout=stdout, stderr=stderr)

        all_charges = list()
        with open(charge_path_main+"/resp2.crg", "r") as all_charges_file:
            for line in all_charges_file:
                l = line.rstrip().split()
                if len(l)==0:
                    continue
                for c in l:
                    all_charges.append(float(c))
        all_charges = np.array(all_charges)
        line_8F10   = fortranformat.FortranRecordWriter('8F10.6')
        atom_counter_old = 0
        atom_counter_new = 0

        ### Finally, write mol2 files for each conformer with the 
        ### calculated resp charges
        for frag_i in range(N_frag):

            mol_frag        = db.get_frag(frag_i)

            if by_charge:
                if charge != Chem.GetFormalCharge(mol_frag):
                    continue

            frag_path_main  = frag_path+"/frag%d" %frag_i

            for conf_i in range(frag_conf_dict[frag_i]):

                charge_path   = frag_path_main+"/conf%d/resp2.crg" %conf_i
                mol2_path     = frag_path_main+"/conf%d/frag%d-conf%d_esp.mol2" \
                                                %(conf_i, frag_i, conf_i)
                mol2_out_path = frag_path_main+"/conf%d/frag%d-conf%d_resp.mol2"\
                                                %(conf_i, frag_i, conf_i)
                mol2_obj      = mol2_reader(mol2_path)
                atom_counter_new = atom_counter_old+mol2_obj.N_atoms

                for i in range(mol2_obj.N_atoms):
                    mol2_obj.charges[i] = all_charges[atom_counter_old:atom_counter_new][i]
                mol2_obj.write(mol2_out_path)

                with open(charge_path, "w") as charge_file:
                    charge_file.write(line_8F10.write(all_charges[atom_counter_old:atom_counter_new]))

                atom_counter_old = atom_counter_new

    ### Finally calculate the charges for the molecules 
    ### that were not fragmented.
    for mol_i, mol_conf in mol_conf_dict.items():

        ### The input files for the resp program will be generated
        ### with the write_resp_in file generator.
        frag_resp_in = write_resp_in()
        esp_list     = list()

        mol_path_main = mol_path+"/frag%d" %mol_i
        mol_path_resp = mol_path_main+"/resp"

        mol_path_esp   = mol_path_resp+"/esp"
        mol_path_resp1 = mol_path_resp+"/resp1"
        mol_path_resp2 = mol_path_resp+"/resp2"

        if not os.path.exists(mol_path_resp):
            os.mkdir(mol_path_resp)

        frag_i_int = 0

        mol_frag   = db.get_mol(mol_i)


        conformer_weights = list()
        if energy_weighting:
            for conf_i in range(mol_conf):
                conf_i_path = mol_path_main+"/conf%d/frag%d-conf%d_esp.log" %(conf_i, mol_i, conf_i)
                energy_list         = list()
                get_energy(conf_i_path, energy_list)
                conformer_weights.append(
                    energy_list[-1] * hartree_to_kcal
                    )
            min_energy = min(conformer_weights)
            for idx in range(len(conformer_weights)):
                conformer_weights[idx] -= min_energy
                ### gas constant in cal/mol/K
                conformer_weights[idx] = np.exp(-conformer_weights[idx]/(1.9872159 * 1.e-3 * 300.))
        else:
            for conf_i in range(mol_conf):
                conformer_weights.append(1.)

        for idx, conf_i in enumerate(range(mol_conf)):
            conf_i_path = mol_path_main+"/conf%d/frag%d-conf%d_esp.log" %(conf_i, mol_i, conf_i)
            esp_path, mol2_path = preprocess_esp(conf_i_path)
            ### Match the mol2 structure obtained from the QM calculation
            ### with the rdkit mol object stored in the database. Perform
            ### some sanity checks on the match lists.
            mol_qm_frag  = Chem.MolFromMol2File(mol2_path, removeHs=False)
            mol_qm_frag  = fix_groups(mol_qm_frag)
            mol_frag     = fix_groups(mol_frag)
            matches      = mol_qm_frag.GetSubstructMatches(mol_frag)
            if len(matches) == 0:
                raise Exception("Did not any substructure matches.")
            if len(matches) > 1:
                raise Exception("Found %d substructure matches. Should be one only." %(len(matches)))
            if len(matches[0]) != mol_frag.GetNumAtoms():
                raise Exception("Found %d atoms in substructure matche. Expected %d." %(len(matches[0]),\
                                 mol_frag.GetNumAtoms()))
            esp_list.append(esp_path)
            frag_resp_in.add_mol(mol_qm_frag,\
                                 Chem.GetFormalCharge(mol_qm_frag),\
                                 conformer_weights[idx],\
                                 "Mol %d Conf %d" %(mol_i,conf_i))
            ### Generate intermolecular restraints
            frag_resp_in.add_intermolecular(0, frag_i_int)

            frag_i_int += 1

        matches = matches[0]

        cat_esp(esp_list, mol_path_esp)
        frag_resp_in.set_title('Stage 1 of multimolecule multiconformation resp fit.')
        frag_resp_in.set_stage('resp1')
        frag_resp_in.write(mol_path_resp1)

        frag_resp_in.set_title('Stage 2 of multimolecule multiconformation resp fit.')
        frag_resp_in.set_stage('resp2')
        frag_resp_in.write(mol_path_resp2)

        resp_args = [resp_exe, "-O", "-i", mol_path_resp1, 
                               "-o", mol_path_resp+"/resp1.out",
                               "-e", mol_path_esp,
                               "-t", mol_path_resp+"/resp1.crg", 
                               "-p", mol_path_resp+"/punch1",
                               "-s", mol_path_resp+"/esout1"]
        print(" ".join(resp_args))
        call(resp_args, stdout=stdout, stderr=stderr)

        resp_args = [resp_exe, "-O", "-i", mol_path_resp2, 
                               "-o", mol_path_resp+"/resp2.out",
                               "-q", mol_path_resp+"/resp1.crg",
                               "-e", mol_path_esp, 
                               "-t", mol_path_resp+"/resp2.crg", 
                               "-p", mol_path_resp+"/punch2",
                               "-s", mol_path_resp+"/esout2"]
        print(" ".join(resp_args))
        call(resp_args, stdout=stdout, stderr=stderr)

        all_charges = list()
        with open(mol_path_resp+"/resp2.crg", "r") as all_charges_file:
            for line in all_charges_file:
                l = line.rstrip().split()
                if len(l)==0:
                    continue
                for c in l:
                    all_charges.append(float(c))
        all_charges = np.array(all_charges)
        line_8F10   = fortranformat.FortranRecordWriter('8F10.6')
        atom_counter_old = 0
        atom_counter_new = 0

        ### Finally, write mol2 files for each conformer with the 
        ### calculated resp charges
        for conf_i in range(mol_conf_dict[mol_i]):

            charge_path   = mol_path_main+"/conf%d/resp2.crg" %conf_i
            mol2_path     = mol_path_main+"/conf%d/frag%d-conf%d_esp.mol2" \
                                            %(conf_i, mol_i, conf_i)
            mol2_out_path = mol_path_main+"/conf%d/frag%d-conf%d_resp.mol2"\
                                            %(conf_i, mol_i, conf_i)
            mol2_obj      = mol2_reader(mol2_path)
            atom_counter_new = atom_counter_old+mol2_obj.N_atoms

            for i in range(mol2_obj.N_atoms):
                mol2_obj.charges[i] = all_charges[atom_counter_old:atom_counter_new][i]
            mol2_obj.write(mol2_out_path)

            with open(charge_path, "w") as charge_file:
                charge_file.write(line_8F10.write(all_charges[atom_counter_old:atom_counter_new]))

            atom_counter_old = atom_counter_new

    ### Map the new charges back to the molecules in the database
    if remap:
        if not os.path.exists(remap_dir):
            os.mkdir(remap_dir)
        for mol_i in range(db.get_mol_count()):
            
            ### First, we must assume, that the Mol object of the molecule
            ### in the database does not include hydrogens. Therefore, we must
            ### match the molecule with full hydrogens on the molecule in our
            ### database.
            mol                         = Chem.MolFromMol2File(db.get_path(mol_i), removeHs=False)
            mol                         = fix_groups(mol)
            db.get_decompose(mol_i).mol = fix_groups(db.get_decompose(mol_i).mol)
            matches_mol                 = mol.GetSubstructMatches(db.get_decompose(mol_i).mol)[0]

            mol_mol2_obj    = mol2_reader(db.get_path(mol_i))
            decomp          = db.get_decompose(mol_i)

            if decomp.get_frag_count() == 0:
                frag_path_main = mol_path+"/frag%d" %mol_i
                frag_mol2_path = frag_path_main+"/conf0/frag%d-conf0_resp.mol2" %mol_i

                mol2_obj_frag  = mol2_reader(frag_mol2_path)

                mol_qm_frag    = Chem.MolFromMol2File(frag_mol2_path, removeHs=False)
                mol_qm_frag    = fix_groups(mol_qm_frag)
                matches_frag   = mol_qm_frag.GetSubstructMatches(mol)[0]

                matched_neighbors = list()
                for atom_idx_frag, atom_idx_mol in zip(matches_frag, matches_mol):

                    mol_mol2_obj.charges[atom_idx_mol] = \
                    mol2_obj_frag.charges[atom_idx_frag]

                    atom_frag = mol_qm_frag.GetAtomWithIdx(atom_idx_frag)
                    atom_mol  = mol.GetAtomWithIdx(atom_idx_mol)
                    for neighbor_frag in atom_frag.GetNeighbors():
                        neighbor_frag_idx = neighbor_frag.GetIdx()
                        if neighbor_frag.GetAtomicNum() == 1:
                            for neighbor_mol in atom_mol.GetNeighbors():
                                neighbor_mol_idx = neighbor_mol.GetIdx()
                                if neighbor_mol.GetAtomicNum() == 1 \
                                and neighbor_mol_idx not in matched_neighbors:
                                    mol_mol2_obj.charges[neighbor_mol_idx] = \
                                    mol2_obj_frag.charges[neighbor_frag_idx]
                                    matched_neighbors.append(neighbor_mol_idx)

            else:
                for frag_i, frag_map, mol_frag in zip(db.get_mol2frag()[mol_i],
                                                      decomp.frag_list_map,
                                                      decomp.frag_list):

                    frag_i_internal = db.get_mol2frag()[mol_i].index(frag_i)

                    frag_path_main  = frag_path+"/frag%d" %frag_i
                    frag_mol2_path  = frag_path_main+"/conf0/frag%d-conf0_resp.mol2" %frag_i

                    mol2_obj_frag   = mol2_reader(frag_mol2_path)

                    ### For each fragment, we must assume, that the atom indexing of the fragment
                    ### as obtained from the decomposition of the molecule is different
                    ### from the indexing of that very fragment in the database.
                    mol_qm_frag     = Chem.MolFromMol2File(frag_mol2_path, removeHs=False)

                    mol_frag        = fix_groups(mol_frag)
                    mol_qm_frag     = fix_groups(mol_qm_frag)
                    matches_frag    = mol_qm_frag.GetSubstructMatches(mol_frag)[0]

                    matched_neighbors = list()

                    for frag_atm_idx, atom_idx_mol in enumerate(frag_map):

                        atom_idx_frag = matches_frag[frag_atm_idx]
                        is_terminal   = False

                        cap_atm_idxs  = list()
                        mol_atm_idxs  = list()

                        for lanc_idxs_count, lanc_idxs in enumerate(decomp.lanc_list[frag_i_internal]):
                            if capsum:
                                break
                            if frag_atm_idx in lanc_idxs:
                                index  = lanc_idxs_count
                                conn_i = decomp.lcap_conn_idx[frag_i_internal][index]
                                conn   = decomp.connectors[conn_i]

                                if not conn.get_terminal():
                                    break
                                for frag_j_internal, rcap_conn_idx in enumerate(decomp.rcap_conn_idx):
                                    if conn_i in rcap_conn_idx:
                                        index_j = rcap_conn_idx.index(conn_i)
                                        ranc_j  = decomp.ranc_list[frag_j_internal][index_j]
                                        for _frag_atm_idx, _atom_idx_mol in enumerate(decomp.frag_list_map[frag_j_internal]):
                                            if _frag_atm_idx in ranc_j:
                                                mol_atm_idxs.append(_atom_idx_mol)

                        for ranc_idxs_count, ranc_idxs in enumerate(decomp.ranc_list[frag_i_internal]):
                            if capsum:
                                break
                            if frag_atm_idx in ranc_idxs:
                                index  = ranc_idxs_count
                                conn_i = decomp.rcap_conn_idx[frag_i_internal][index]
                                conn   = decomp.connectors[conn_i]

                                if not conn.get_terminal():
                                    break

                                for frag_j_internal, lcap_conn_idx in enumerate(decomp.lcap_conn_idx):
                                    if conn_i in lcap_conn_idx:
                                        index_j = lcap_conn_idx.index(conn_i)
                                        lanc_j  = decomp.lanc_list[frag_j_internal][index_j]
                                        for _frag_atm_idx, _atom_idx_mol in enumerate(decomp.frag_list_map[frag_j_internal]):
                                            if _frag_atm_idx in lanc_j:
                                                #print (frag_i_internal,)
                                                #print (index_j,)
                                                #print (lanc_j,)
                                                #print (db.get_name(mol_i), frag_i,)
                                                #print (db.get_mol2frag()[mol_i][frag_i_internal],)
                                                #print (frag_j_internal,)
                                                #print (_atom_idx_mol)
                                                mol_atm_idxs.append(_atom_idx_mol)

                        if is_terminal and not capsum:
                            if len(cap_atm_idxs) != len(mol_atm_idxs):
                                print (cap_atm_idxs)
                                print (mol_atm_idxs)
                                raise Exception("Could not find all corresponding cap atoms in molecule during non-capsum mode. \
Either check your linker definition or run in capsum mode.")

                        for i in range(len(cap_atm_idxs)):
                            ### initial rdkit object fragment -> fragment after QM calculation with hydrogen
                            cap_atm_idxs[i] = matches_frag[cap_atm_idxs[i]]
                            mol_atm_idxs[i] = matches_mol[mol_atm_idxs[i]]

                        #print (db.get_name(mol_i), frag_i,)
                        #print (db.get_mol2frag()[mol_i][frag_i_internal],)
                        #print (cap_atm_idxs,)
                        #print (mol_atm_idxs)

                        ### Match charges from fragment to molecule for atoms
                        ### that were already present in the initial molecule.
                        ### Note that we must assume, that hydrogens are likely
                        ### not present in the initial molecule.
                        mol_mol2_obj.charges[matches_mol[atom_idx_mol]] = \
                        mol2_obj_frag.charges[atom_idx_frag]

                        atom_frag = mol_qm_frag.GetAtomWithIdx(atom_idx_frag)
                        atom_mol  = mol.GetAtomWithIdx(matches_mol[atom_idx_mol])
                        for neighbor_frag in atom_frag.GetNeighbors():
                            neighbor_frag_idx = neighbor_frag.GetIdx()
                            if neighbor_frag.GetAtomicNum() == 1 \
                            and neighbor_frag_idx not in cap_atm_idxs:
                                for neighbor_mol in atom_mol.GetNeighbors():
                                    neighbor_mol_idx = neighbor_mol.GetIdx()
                                    if neighbor_mol.GetAtomicNum() == 1 \
                                    and neighbor_mol_idx not in matched_neighbors:
                                        mol_mol2_obj.charges[neighbor_mol_idx] = \
                                        mol2_obj_frag.charges[neighbor_frag_idx]
                                        matched_neighbors.append(neighbor_mol_idx)

                            elif is_terminal and not capsum \
                            and neighbor_frag_idx in cap_atm_idxs:
                                neighbor_mol_idx = mol_atm_idxs[cap_atm_idxs.index(neighbor_frag_idx)]
                                if neighbor_mol_idx in matched_neighbors:
                                    continue
                                mol_mol2_obj.charges[neighbor_mol_idx] = \
                                mol2_obj_frag.charges[neighbor_frag_idx]
                                matched_neighbors.append(neighbor_mol_idx)

                                for neighbor_frag2 in neighbor_frag.GetNeighbors():
                                    neighbor_frag_idx2 = neighbor_frag2.GetIdx()
                                    if neighbor_frag2.GetAtomicNum() == 1:
                                        atom_mol2 = mol.GetAtomWithIdx(neighbor_mol_idx)
                                        for neighbor_mol2 in atom_mol2.GetNeighbors():
                                            neighbor_mol_idx2 = neighbor_mol2.GetIdx()
                                            if neighbor_mol2.GetAtomicNum() == 1 \
                                            and neighbor_mol_idx2 not in matched_neighbors:
                                                mol_mol2_obj.charges[neighbor_mol_idx2] = \
                                                mol2_obj_frag.charges[neighbor_frag_idx2]
                                                matched_neighbors.append(neighbor_mol_idx2)

                                                #print (mol_mol2_obj.charges[neighbor_mol_idx], )
                                                #print (neighbor_mol_idx)
                                                #print ("->", db.get_name(mol_i), )
                                                #print (mol.GetAtomWithIdx(neighbor_mol_idx).GetAtomicNum(), )
                                                #print (frag_i, )
                                                #print (neighbor_mol_idx, )
                                                #print (neighbor_frag_idx)

                        #print (mol_mol2_obj.charges[12])

                #print ("Finished", db.get_name(mol_i))
                mol_mol2_obj.write(remap_dir+"/"+db.get_name(mol_i)+".resp.mol2")