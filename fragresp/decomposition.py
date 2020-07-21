import copy
from rdkit import Chem

from fragresp.utils import are_conn_equal
from fragresp.utils import get_parity
from fragresp.utils import get_frag_bonds
from fragresp.utils import canonicalize_tautomers
from fragresp.constants import CHI_TETRAHEDRAL_CW
from fragresp.constants import CHI_TETRAHEDRAL_CCW

class decompose(object):

    def __init__(self, Mol, Verbose=False):

        self.mol                = Mol

        self.canonical_rank     = list()
        string_rank             = list(Chem.CanonicalRankAtoms(self.mol, breakTies=False))
        for rank in string_rank:
            self.canonical_rank.append(int(rank))
        del string_rank
        canonicalize_tautomers(self.canonical_rank, self.mol)

        ### Holds connector instances
        self.connectors         = list()
        
        ### Holds rdkit Mol instancs of final capped fragments
        self.frag_list          = list()
        ### Holds atom indices of fragments in numbering scheme of
        ### original molecule.
        self.frag_list_map      = list()

        ### Holds r/l anchor atom idcs for each fragment in numbering
        ### scheme of the fragment molecule
        self.ranc_list          = list()
        self.lanc_list          = list()

        ### Holds r/l cap atom idcs for each fragment in numbering
        ### scheme of the fragment molecule
        self.rcap_list_map      = list()
        self.lcap_list_map      = list()

        ### Holds corresponding connector idx for each r/l cap
        self.rcap_conn_idx      = list()
        self.lcap_conn_idx      = list()

        ### Stores fragment to fragment cross couplings
        ### atom indices for each cross couplin
        self.frag2frag_atms  = list()
        ### fragment indices for each cross coulin
        self.frag2frag_frgs  = list()

        self.__frag_count       = 0

        self.__connector_count  = 0
        
        self.verbose            = Verbose
        
        if self.verbose:
            self.process_list = list()

    def get_frag_count(self):
        return self.__frag_count

    def add_connector(self, connector):

        #, lcap_smi, rcap_smi, lanc=None, ranc=None, ring=None):

        if self.verbose:
            print ("Adding connctor", self.__connector_count)

        if connector.check():

            if self.__connector_count > 0:

                for conn in self.connectors:

                    if self.verbose:
                        print ("Comparing", connector.get_name(), "with", conn.get_name())

                    if are_conn_equal(connector, conn):
                        raise Warning ("Connector already present.")
                        return False

            self.connectors.append(connector)

            self.__connector_count += 1

            return True

        else:

            raise Warning("Connector not fully defined.")
            return False

    def decompose(self):

        ### Carry out the actual fragmentization
        if self.verbose:
            print ("Start molecule decomposition...")

        if self.__decompose_mol(self.mol)==True:
            if self.verbose:
                print ("Done.")
                print ("")
            return True
        else:
            return False

    def __decompose_mol(self, Mol):

        ### Mol: Molecule that will be decompsed. Must be an instance of
        ###      Rdkit molecule object.

        ### Holds bond indices of broken bonds at r/l anchor
        ranc_bond_idcs = list()
        lanc_bond_idcs = list()

        ### Holds atom indices of broken bonds at r/l anchor
        ranc_atom_idcs = list()
        lanc_atom_idcs = list()

        ### Holds atom indices of bonds at terminal r/l anchor. 
        ### They won't be broken.
        ranc_atom_idcs_t = list()
        lanc_atom_idcs_t = list()

        ### Holds atom indices of broken bonds at connectors
        rconn_atom_idcs = list()
        lconn_atom_idcs = list()

        ### Holds atom indices of bonds between fragment and
        ### terminal connectors. They won't be broken.
        rconn_atom_idcs_t = list()
        lconn_atom_idcs_t = list()

        ### Holds permutation parity of atoms involved in broken bonds
        ranc_parity = list()
        lanc_parity = list()

        ### Stores if anchor atom is chiral or not
        ranc_chiral = list()
        lanc_chiral = list()

        ### Stores bond type information for bonds between fragment and connector/cap
        ranc_bondtype  = list()
        lanc_bondtype  = list()

        ### Stores the connector index of the bond
        rconn_idcs = list()
        lconn_idcs = list()

        ### Stores the connector index of the bond with terminal anchor
        rconn_idcs_t = list()
        lconn_idcs_t = list()

        ### Stores atom indices of all anchor atoms in molecule
        ranc_mol = list()
        lanc_mol = list()

        ### Stores atom indices of all terminal anchor atoms in molecule 
        ranc_mol_t = list()
        lanc_mol_t = list()

        ### List that holds all atom indices
        Mol_atm_idxs = range(Mol.GetNumAtoms())

        if self.verbose:
            print ("Start analyzing bonds...")

        ### Begin loop over all connectors
        for connector_idx in range(self.__connector_count):

            conn_obj  = self.connectors[connector_idx]

            connector = conn_obj.get_connector()
            lanc      = conn_obj.get_lanc()
            ranc      = conn_obj.get_ranc()
            ring      = conn_obj.get_ring()

            ranc_map  = conn_obj.get_ranc_map()
            lanc_map  = conn_obj.get_lanc_map()

            connector_mol_matches       = Mol.GetSubstructMatches(connector, useChirality=True)
            connector_mol_match_count   = len(connector_mol_matches)

            connector_mol_match_exclude = list()

            if self.verbose:
                print ("Decomposing with connector ID %d (%s)..." %(connector_idx, self.connectors[connector_idx].get_name()))
                print ("Found %d matches." %connector_mol_match_count)

            if connector_mol_match_count == 0:
                if self.verbose:
                    print ("No matches found for current connector.",)
                    print ("No decomposition.")
                continue

            ### Find the corresponding atom indices of the anchor atoms
            lanc_mol_matches = list()
            ranc_mol_matches = list()
            for idx1 in range(connector_mol_match_count):
                ranc_mol_matches.append(list())
                lanc_mol_matches.append(list())
                for r_i in ranc:
                    ranc_mol_matches[-1].append(connector_mol_matches[idx1][r_i])
                for l_i in lanc:
                    lanc_mol_matches[-1].append(connector_mol_matches[idx1][l_i])

            ### Filter substructure matches that do not unambigously
            ### identify connectors, i.e. they have the same connector 
            ### atoms except the anchor atoms. Also, identify all
            ### connectors, that overlap with other, previously matched,
            ### connectors.
            ### Add those connector matches to an exclcude list.
            for idx1 in range(connector_mol_match_count):
                for idx2 in range(connector_mol_match_count):
                    if idx1<idx2:
                        for idx1_match in connector_mol_matches[idx1]:
                            idx1_match_1 = -1
                            idx1_match_2 = -1
                            if idx1_match in lanc_mol_matches[idx1]:
                                idx1_match_1 = lanc_mol_matches[idx1].index(idx1_match)
                            if idx1_match in lanc_mol_matches[idx2]:
                                idx1_match_2 = lanc_mol_matches[idx2].index(idx1_match)
                            if idx1_match_1>-1 and idx1_match_1==idx1_match_2:
                                if lanc_map[idx1_match_1]==-1:
                                    continue
                                if idx1 not in connector_mol_match_exclude:
                                    connector_mol_match_exclude.append(idx1)
                                if idx2 not in connector_mol_match_exclude:
                                    connector_mol_match_exclude.append(idx2)
                                if self.verbose:
                                    print ("Excluded connector match due to",)
                                    print ("ambigious matching on lanc atoms",)
                                    print ("%s (match id %d,%d)." %(idx1_match, idx1, idx2))
                                continue

                            idx1_match_1 = -1
                            idx1_match_2 = -1
                            if idx1_match in ranc_mol_matches[idx1]:
                                idx1_match_1 = ranc_mol_matches[idx1].index(idx1_match)
                            if idx1_match in ranc_mol_matches[idx2]:
                                idx1_match_2 = ranc_mol_matches[idx2].index(idx1_match)
                            if idx1_match_1>-1 and idx1_match_1==idx1_match_2:
                                if ranc_map[idx1_match_1]==-1:
                                    continue
                                if idx1 not in connector_mol_match_exclude:
                                    connector_mol_match_exclude.append(idx1)
                                if idx2 not in connector_mol_match_exclude:
                                    connector_mol_match_exclude.append(idx2)
                                if self.verbose:
                                    print ("Excluded connector match due to",)
                                    print ("ambigious matching on ranc atoms",)
                                    print ("%s (match id %d,%d)." %(idx1_match, idx1, idx2))
                                continue

            if self.verbose:
                print ("Excluded connector matches: ", connector_mol_match_exclude)

            ### Begin loop over connector_mol_matches
            for idx1 in range(connector_mol_match_count):
                connector_mol_match = connector_mol_matches[idx1]
                ranc_mol_match      = ranc_mol_matches[idx1]
                lanc_mol_match      = lanc_mol_matches[idx1]

                if conn_obj.get_terminal():
                    for term_anc_idx in conn_obj.get_terminal_anc():
                        term_anc_atm = Mol.GetAtomWithIdx(connector_mol_match[term_anc_idx])
                        for term_neighbor in term_anc_atm.GetNeighbors():
                            term_neighbor_idx = term_neighbor.GetIdx()
                            if term_neighbor_idx not in connector_mol_match:
                                connector_mol_match_exclude.append(idx1)
                                if self.verbose:
                                    print ("Excluded current connector match",)
                                    print ("due to terminal anchor atom %d at non-terminal" %connector_mol_match[term_anc_idx],)
                                    print ("position connected to atom %d." %term_neighbor_idx)

                if self.verbose:
                    print ("Current connector_mol_match ID", idx1)
                    print ("Current connector:", connector_mol_match)

                ranc_bond_idcs_tmp, ratom_idcs_tmp = get_frag_bonds(Mol,
                                                                    ranc_mol_match, 
                                                                    connector_mol_match,
                                                                    ring,
                                                                    self.verbose)

                lanc_bond_idcs_tmp, latom_idcs_tmp = get_frag_bonds(Mol, 
                                                                    lanc_mol_match, 
                                                                    connector_mol_match, 
                                                                    ring,
                                                                    self.verbose)

                _check_ranc_atom_idcs  = ranc_atom_idcs  + ranc_atom_idcs_t
                _check_lanc_atom_idcs  = lanc_atom_idcs  + lanc_atom_idcs_t
                _check_rconn_atom_idcs = rconn_atom_idcs + rconn_atom_idcs_t
                _check_lconn_atom_idcs = lconn_atom_idcs + lconn_atom_idcs_t
                _check_rconn_idcs      = rconn_idcs      + rconn_idcs_t
                _check_lconn_idcs      = lconn_idcs      + lconn_idcs_t
                _check_ranc_mol        = ranc_mol        + ranc_mol_t
                _check_lanc_mol        = lanc_mol        + lanc_mol_t
                ### Sanity checking for ranc
                for _anc_atom, _anc_bond in zip(ratom_idcs_tmp,ranc_bond_idcs_tmp):
                    for list_idx in range(len(_check_ranc_mol)):
                        _ranc_mol_match = _check_ranc_mol[list_idx]
                        _lanc_mol_match = _check_lanc_mol[list_idx]
                        rconnector_idx  = _check_rconn_idcs[list_idx]
                        lconnector_idx  = _check_lconn_idcs[list_idx]

                        if _check_ranc_atom_idcs[list_idx] in ranc_mol_match \
                        and _check_rconn_atom_idcs[list_idx] in ranc_mol_match:
                            if self.verbose:
                                print ("Excluded current connector match",)
                                print ("due to ranc overlap with ranc of connector", rconnector_idx,)
                                print (self.connectors[rconnector_idx].get_name())
                            if idx1 not in connector_mol_match_exclude:
                                connector_mol_match_exclude.append(idx1)

                        elif _check_lanc_atom_idcs[list_idx] in ranc_mol_match \
                        and _check_lconn_atom_idcs[list_idx] in ranc_mol_match:
                            if self.verbose:
                                print ("Excluded current connector match",)
                                print ("due to ranc overlap with lanc of connector", lconnector_idx,)
                                print (self.connectors[lconnector_idx].get_name())
                            if idx1 not in connector_mol_match_exclude:
                                connector_mol_match_exclude.append(idx1)

                        elif (_check_ranc_atom_idcs[list_idx]==_anc_atom[0] \
                        and _check_rconn_atom_idcs[list_idx]==_anc_atom[1]) \
                        or (_check_ranc_atom_idcs[list_idx]==_anc_atom[1] \
                        and _check_rconn_atom_idcs[list_idx]==_anc_atom[0]):
                            if self.verbose:
                                print ("Excluded current connector match",)
                                print ("due to ranc overlap with ranc of connector", rconnector_idx,)
                                print (self.connectors[rconnector_idx].get_name())
                            if idx1 not in connector_mol_match_exclude:
                                connector_mol_match_exclude.append(idx1)

                        elif (_check_lanc_atom_idcs[list_idx]==_anc_atom[0] \
                        and _check_lconn_atom_idcs[list_idx]==_anc_atom[1]) \
                        or (_check_lanc_atom_idcs[list_idx]==_anc_atom[1] \
                        and _check_lconn_atom_idcs[list_idx]==_anc_atom[0]):
                            if self.verbose:
                                print ("Excluded current connector match",)
                                print ("due to ranc overlap with lanc of connector", lconnector_idx,)
                                print (self.connectors[lconnector_idx].get_name())
                            if idx1 not in connector_mol_match_exclude:
                                connector_mol_match_exclude.append(idx1)

                        elif _anc_atom[0] in _ranc_mol_match \
                        and _anc_atom[1] in _ranc_mol_match:
                            if self.verbose:
                                print ("Excluded current connector match",)
                                print ("due to ranc overlap with ranc of connector", rconnector_idx,)
                                print (self.connectors[rconnector_idx].get_name())
                            if idx1 not in connector_mol_match_exclude:
                                connector_mol_match_exclude.append(idx1)

                        elif _anc_atom[0] in _lanc_mol_match \
                        and _anc_atom[1] in _lanc_mol_match:
                            if self.verbose:
                                print ("Excluded current connector match",)
                                print ("due to ranc overlap with lanc of connector", lconnector_idx,)
                                print (self.connectors[lconnector_idx].get_name())
                            if idx1 not in connector_mol_match_exclude:
                                connector_mol_match_exclude.append(idx1)

                    if self.verbose and idx1 not in connector_mol_match_exclude:
                        print ("Ranc decompostion of bond %d of anc atom %d - conn atom %d." %(_anc_bond,
                                                                                            _anc_atom[0], 
                                                                                            _anc_atom[1]))

                ### Sanity checking for lanc
                for _anc_atom, _anc_bond in zip(latom_idcs_tmp,lanc_bond_idcs_tmp):
                    for list_idx in range(len(lanc_mol)):
                        _ranc_mol_match = _check_ranc_mol[list_idx]
                        _lanc_mol_match = _check_lanc_mol[list_idx]
                        rconnector_idx  = _check_rconn_idcs[list_idx]
                        lconnector_idx  = _check_lconn_idcs[list_idx]

                        if _check_ranc_atom_idcs[list_idx] in ranc_mol_match \
                        and _check_rconn_atom_idcs[list_idx] in ranc_mol_match:
                            if self.verbose:
                                print ("Excluded current connector match",)
                                print ("due to lanc overlap with ranc of connector", rconnector_idx,)
                                print (self.connectors[rconnector_idx].get_name())
                            if idx1 not in connector_mol_match_exclude:
                                connector_mol_match_exclude.append(idx1)

                        elif _check_lanc_atom_idcs[list_idx] in ranc_mol_match \
                        and _check_lconn_atom_idcs[list_idx] in ranc_mol_match:
                            if self.verbose:
                                print ("Excluded current connector match",)
                                print ("due to lanc overlap with lanc of connector", lconnector_idx,)
                                print (self.connectors[lconnector_idx].get_name())
                            if idx1 not in connector_mol_match_exclude:
                                connector_mol_match_exclude.append(idx1)

                        elif (_check_ranc_atom_idcs[list_idx]==_anc_atom[0] \
                        and _check_rconn_atom_idcs[list_idx]==_anc_atom[1]) \
                        or (_check_ranc_atom_idcs[list_idx]==_anc_atom[1] \
                        and _check_rconn_atom_idcs[list_idx]==_anc_atom[0]):
                            if self.verbose:
                                print ("Excluded current connector match",)
                                print ("due to lanc overlap with ranc of connector", rconnector_idx,)
                                print (self.connectors[rconnector_idx].get_name())
                            if idx1 not in connector_mol_match_exclude:
                                connector_mol_match_exclude.append(idx1)

                        elif (_check_lanc_atom_idcs[list_idx]==_anc_atom[0] \
                        and _check_lconn_atom_idcs[list_idx]==_anc_atom[1]) \
                        or (_check_lanc_atom_idcs[list_idx]==_anc_atom[1] \
                        and _check_lconn_atom_idcs[list_idx]==_anc_atom[0]):
                            if self.verbose:
                                print ("Excluded current connector match",)
                                print ("due to lanc overlap with lanc of connector", lconnector_idx,)
                                print (self.connectors[lconnector_idx].get_name())
                            if idx1 not in connector_mol_match_exclude:
                                connector_mol_match_exclude.append(idx1)

                        elif _anc_atom[0] in _ranc_mol_match \
                        and _anc_atom[1] in _ranc_mol_match:
                            if self.verbose:
                                print ("Excluded current connector match",)
                                print ("due to lanc overlap with ranc of connector", rconnector_idx,)
                                print (self.connectors[rconnector_idx].get_name())
                            if idx1 not in connector_mol_match_exclude:
                                connector_mol_match_exclude.append(idx1)

                        elif _anc_atom[0] in _lanc_mol_match \
                        and _anc_atom[1] in _lanc_mol_match:
                            if self.verbose:
                                print ("Excluded current connector match",)
                                print ("due to lanc overlap with lanc of connector", lconnector_idx,)
                                print (self.connectors[lconnector_idx].get_name())
                            if idx1 not in connector_mol_match_exclude:
                                connector_mol_match_exclude.append(idx1)

                    if self.verbose and idx1 not in connector_mol_match_exclude:
                        print ("Lanc decompostion of bond %d of anc atom %d - conn atom %d." %(_anc_bond,
                                                                                            _anc_atom[0], 
                                                                                            _anc_atom[1]))

                if len(ranc_bond_idcs_tmp)==0 and \
                   len(lanc_bond_idcs_tmp)==0:
                    connector_mol_match_exclude.append(idx1)

                if idx1 in connector_mol_match_exclude:
                    if self.verbose:
                        print ("Current connector_mol_match excluded.")
                    continue

                if self.verbose:
                    print ("R anchor atom indices in molecule:", ranc_mol_match)
                    print ("L anchor atom indices in molecule:", lanc_mol_match)

                ### Keep track of bonds that will be broken (r site)
                for bond_idx,\
                    atom_idx in zip(ranc_bond_idcs_tmp,
                                    ratom_idcs_tmp):

                    anc_idx, conn_idx = atom_idx

                    atm_anc  = Mol.GetAtomWithIdx(anc_idx)
                    atm_conn = Mol.GetAtomWithIdx(conn_idx)

                    if conn_obj.get_terminal():
                        ranc_atom_idcs_t.append(anc_idx)
                        ranc_mol_t.append(ranc_mol_match)
                        rconn_atom_idcs_t.append(conn_idx)
                        rconn_idcs_t.append(connector_idx)

                    else:
                        ranc_bond_idcs.append(bond_idx)

                        ranc_atom_idcs.append(anc_idx)
                        ranc_mol.append(ranc_mol_match)
                        rconn_atom_idcs.append(conn_idx)
                        rconn_idcs.append(connector_idx)

                ### Keep track of bonds that will be broken (l site)
                for bond_idx,\
                    atom_idx in zip(lanc_bond_idcs_tmp,
                                    latom_idcs_tmp):

                    anc_idx, conn_idx = atom_idx

                    atm_anc  = Mol.GetAtomWithIdx(anc_idx)
                    atm_conn = Mol.GetAtomWithIdx(conn_idx)

                    if conn_obj.get_terminal():
                        lanc_atom_idcs_t.append(anc_idx)
                        lanc_mol_t.append(lanc_mol_match)
                        lconn_atom_idcs_t.append(conn_idx)
                        lconn_idcs_t.append(connector_idx)

                    else:
                        lanc_bond_idcs.append(bond_idx)

                        lanc_atom_idcs.append(anc_idx)
                        lanc_mol.append(lanc_mol_match)
                        lconn_atom_idcs.append(conn_idx)
                        lconn_idcs.append(connector_idx)

            ### End loop over connector_mol_matches

        ### End loop over all connectors

        if len(ranc_bond_idcs)>0 or len(lanc_bond_idcs)>0:

            if self.verbose:
                print ("Start decomposing molecule....")
            
            frags = Chem.FragmentOnBonds(Mol,
                                         list(set(ranc_bond_idcs+lanc_bond_idcs)), #Must be uniq
                                         addDummies=False)

            frags_mol, mol_frags_atm_idxs = Chem.GetMolFrags(frags, asMols=True, sanitizeFrags=False),\
                                            Chem.GetMolFrags(frags, asMols=False, sanitizeFrags=False)

            if self.verbose:
                self.process_list = frags_mol

            ### Make fragment to fragment cross couplings
            frag_canonical_ranks = list()
            for frag_mol, mol_atm_idxs in zip(frags_mol, mol_frags_atm_idxs):
                frag_canonical_ranks.append(list())
                for atm_idx, mol_atm_idx in enumerate(mol_atm_idxs):
                    frag_canonical_ranks[-1].append(self.canonical_rank[mol_atm_idx])

            N_frags = len(frag_canonical_ranks)
            for i1 in range(N_frags):
                canonical_ranks1 = frag_canonical_ranks[i1]
                for i2 in range(N_frags):
                    if (i2-1)<i1:
                        continue
                    canonical_ranks2 = frag_canonical_ranks[i2]
                    for atm_idx2, rank2 in enumerate(canonical_ranks2):
                        if rank2 in canonical_ranks1:
                            atm_idx1 = canonical_ranks1.index(rank2)
                            self.frag2frag_frgs.append([i1,i2])
                            self.frag2frag_atms.append([atm_idx1,atm_idx2])

            ###
            ### Loop over all fragments generated by FragmentOnBonds.
            ### For each fragment that is not a connector find the 
            ### corresponding connector fragment.
            ###
            ### frag_mol    : rdkit Mol instance
            ### mol_atm_idxs: atom indices of the original molecule that
            ###               generated the Mol instance of the original molecule
            ###

            for frag_mol, mol_atm_idxs in zip(frags_mol, mol_frags_atm_idxs):

                atm_idx_list   = list()
                parity_list    = list()
                chirality_list = list()
                bondtype_list  = list()

                combo          = copy.copy(frag_mol)

                rcap_list_map = list()
                lcap_list_map = list()

                rcap_conn_idx = list()
                lcap_conn_idx = list()

                ranc_list     = list()
                lanc_list     = list()

                ###
                ### Loop over all atoms in the fragment
                ###
                ### atm_idx    : atom index in the fragment indexing scheme
                ### mol_atm_idx: atom index in the original molecule indexing scheme
                ###
                ### Find fragments and the corresponding connector entry
                ###
                for atm_idx, mol_atm_idx in enumerate(mol_atm_idxs):

                    ###
                    ### A connector atom is stored in connector_mol_idcs
                    ### but will not be found in the anchor lists ranc_atom_idcs 
                    ### and lanc_atom_idcs.
                    ###

                    ### Work through the terminal groups.
                    ### They have not been modified, since no
                    ### bonds are broken here.
                    for list_idx, _ranc_idx in enumerate(ranc_atom_idcs_t):
                        if mol_atm_idx != _ranc_idx:
                            continue

                        connector_idx = rconn_idcs_t[list_idx]

                        ranc_list.append(list())
                        for _anc_idx in ranc_mol_t[list_idx]:
                            if _anc_idx in mol_atm_idxs:
                                ranc_list[-1].append(mol_atm_idxs.index(_anc_idx))

                        rcap_conn_idx.append(connector_idx)
                        rcap_list_map.append(list())

                    for list_idx, _lanc_idx in enumerate(lanc_atom_idcs_t):
                        if mol_atm_idx != _lanc_idx:
                            continue

                        connector_idx = lconn_idcs_t[list_idx]

                        lanc_list.append(list())
                        for _anc_idx in lanc_mol_t[list_idx]:
                            if _anc_idx in mol_atm_idxs:
                                lanc_list[-1].append(mol_atm_idxs.index(_anc_idx))

                        lcap_conn_idx.append(connector_idx)
                        lcap_list_map.append(list())

                    ### Now, work through the non-terminal groups.
                    ### All these groups have open bond endings, since we have
                    ### broken them during the fragmentation process. Connect
                    ### all these open bond endings to capping groups.
                    for list_idx, _ranc_idx in enumerate(ranc_atom_idcs):
                        if mol_atm_idx != _ranc_idx:
                            continue

                        ranc_list.append(list())
                        for _anc_idx in ranc_mol[list_idx]:
                            if _anc_idx in mol_atm_idxs:
                                ranc_list[-1].append(mol_atm_idxs.index(_anc_idx))

                        ### combo_length is also the offset of the idx for
                        ### the newly added cap.
                        connector_idx = rconn_idcs[list_idx]
                        conn          = self.connectors[connector_idx]
                        cap           = conn.get_rcap()
                        combo_length  = combo.GetNumAtoms()
                        cap_length    = cap.GetNumAtoms()

                        ### Connect fragment and connector cap
                        combo = Chem.CombineMols(combo, cap)
                        atm_idx_list.append([atm_idx,
                                             combo_length])
                        parity_list.append([get_parity(Mol, mol_atm_idx), 
                                            get_parity(cap, 0)])
                        chirality_list.append([Mol.GetAtomWithIdx(atm_idx).GetChiralTag(),
                                               cap.GetAtomWithIdx(0).GetChiralTag()])
                        bondtype_list.append(Mol.GetBondWithIdx(ranc_bond_idcs[list_idx]).GetBondType())

                        rcap_list_map.append(range(combo_length, combo_length+cap_length))
                        rcap_conn_idx.append(connector_idx)


                    for list_idx, _lanc_idx in enumerate(lanc_atom_idcs):

                        if mol_atm_idx != _lanc_idx:
                            continue

                        lanc_list.append(list())
                        for _anc_idx in lanc_mol[list_idx]:
                            if _anc_idx in mol_atm_idxs:
                                lanc_list[-1].append(mol_atm_idxs.index(_anc_idx))

                        ### combo_length is also the offset of the idx for
                        ### the newly added cap.
                        connector_idx = lconn_idcs[list_idx]
                        conn          = self.connectors[connector_idx]
                        cap           = conn.get_lcap()
                        combo_length  = combo.GetNumAtoms()
                        cap_length    = cap.GetNumAtoms()

                        ### Connect fragment and connector cap
                        combo = Chem.CombineMols(combo, cap)
                        atm_idx_list.append([atm_idx,
                                             combo_length])
                        parity_list.append([get_parity(Mol, mol_atm_idx), 
                                            get_parity(cap, 0)])
                        chirality_list.append([Mol.GetAtomWithIdx(atm_idx).GetChiralTag(),
                                               cap.GetAtomWithIdx(0).GetChiralTag()])
                        bondtype_list.append(Mol.GetBondWithIdx(lanc_bond_idcs[list_idx]).GetBondType())

                        lcap_list_map.append(range(combo_length, combo_length+cap_length))
                        lcap_conn_idx.append(connector_idx)

                ### Put molecule together
                e_combo = Chem.RWMol(combo)
                for atm_idx, bondtype in zip(atm_idx_list, bondtype_list):
                    e_combo.AddBond(atm_idx[0], atm_idx[1], bondtype)

                combo = e_combo.GetMol()

                for atm_idx, parity, chirality in zip(atm_idx_list,\
                                                      parity_list,\
                                                      chirality_list):

                    for i in range(2):
                        if get_parity(combo, atm_idx[i]) != parity[i]:
                            new_chiral = None
                            if chirality[i] == CHI_TETRAHEDRAL_CCW:
                                new_chiral = CHI_TETRAHEDRAL_CW
                            elif chirality[i] == CHI_TETRAHEDRAL_CW:
                                new_chiral = CHI_TETRAHEDRAL_CCW
                            if new_chiral != None:
                                atm = e_combo.GetAtomWithIdx(atm_idx[i])
                                atm.SetChiralTag(new_chiral)

                self.frag_list.append(combo)
                self.frag_list_map.append(mol_atm_idxs)

                self.rcap_list_map.append(rcap_list_map)
                self.lcap_list_map.append(lcap_list_map)

                self.rcap_conn_idx.append(rcap_conn_idx)
                self.lcap_conn_idx.append(lcap_conn_idx)

                self.ranc_list.append(ranc_list)
                self.lanc_list.append(lanc_list)

                self.__frag_count += 1

            ### end loop over all fragments

        elif self.verbose:
            print ("Nothing to decompose.")

        return True
