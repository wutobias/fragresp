from rdkit import Chem
import copy

from fragresp.utils import are_mol_same, are_conn_equal, get_parity
from fragresp.constants import CHI_TETRAHEDRAL_CW
from fragresp.constants import CHI_TETRAHEDRAL_CCW

class fragment_db(object):
    
    def __init__(self, Verbose=False):

        self.__mol_path  = list()
        self.__name_list = list()

        self.__dec_list      = list()
        self.__mol2frag      = list()
        self.__mol2fragmap   = list()

        self.__frag_list     = list()
        self.__frag2mol      = list()

        self.__conn_list     = list()
        self.__conn2mol      = list()
        self.__mol2conn      = list()

        self.__mol_count     = 0
        self.__frag_count    = 0
        self.__conn_count    = 0

        self.__frag2rcapconn = list()
        self.__frag2lcapconn = list()

        self.__frag2rcapmap  = list()
        self.__frag2lcapmap  = list()

        self.__frag2rancmap  = list()
        self.__frag2lancmap  = list()

        self.__frag2rancmap_t = list()
        self.__frag2lancmap_t = list()

        self.__frag2rancconn_t = list()
        self.__frag2lancconn_t = list()

        self.__frag2frag_frgs = list()
        self.__frag2frag_atms = list()

        self.verbose = Verbose
        
    def __incr_mol(self):
        self.__mol_count += 1
        
    def __incr_frag(self):
        self.__frag_count += 1

    def __incr_conn(self):
        self.__conn_count += 1

    def get_mol_count(self):
        return self.__mol_count

    def get_frag_count(self):
        return self.__frag_count

    def get_conn_count(self):
        return self.__conn_count

    def get_path_list(self):
        return self.__mol_path
    
    def get_path(self, i):
        return self.__mol_path[i]

    def get_name_list(self):
        return self.__name_list
    
    def get_name(self, i):
        return self.__name_list[i]
    
    def get_frag_list(self):
        return self.__frag_list
    
    def get_frag(self, i):
        return self.__frag_list[i]

    def get_mol2frag(self):
        return self.__mol2frag

    def get_frag2mol(self):
        return self.__frag2mol

    def get_mol(self,i):
        return self.__dec_list[i].mol

    def get_mol_list(self):
        return [s.mol for s in self.__dec_list]

    def get_surr(self, i):
        if not self.__conn_list[i].get_terminal():
            return self.__conn_list[i].get_surrogate()
        return None

    def get_surr_list(self):
        surr_list = list()
        for i in range(self.__conn_count):
            surr_list.append(self.get_surr(i))
        return surr_list

    def get_surr_cap(self, i):
        if not self.__conn_list[i].get_terminal():
            return self.__conn_list[i].get_surrogate_cap()
        return None

    def get_surr_cap_list(self):
        surr_cap_list = list()
        for i in range(self.__conn_count):
            surr_cap_list.append(self.get_surr_cap(i))
        return surr_cap_list

    def get_decompose_list(self):
        return self.__dec_list

    def get_decompose(self, i):
        return self.__dec_list[i]

    def get_mol2fragmap(self):
        return self.__mol2fragmap

    def get_mol2conn(self):
        return self.__mol2conn

    def get_conn2mol(self):
        return self.__conn2mol

    def get_conn(self,i):
        return self.__conn_list[i]

    def get_conn_list(self):
        return self.__conn_list

    def get_frag2rcapconn(self):
        return self.__frag2rcapconn

    def get_frag2lcapconn(self):
        return self.__frag2lcapconn

    def get_frag2lcapmap(self):
        return self.__frag2lcapmap

    def get_frag2rcapmap(self):
        return self.__frag2rcapmap

    def get_frag2rancmap(self):
        return self.__frag2rancmap

    def get_frag2lancmap(self):
        return self.__frag2lancmap

    def get_frag2frag_frgs(self):
        return self.__frag2frag_frgs

    def get_frag2frag_atms(self):
        return self.__frag2frag_atms

    def add_mol(self,
                decompose,
                path,
                name):

        if self.verbose:
            print ("Current parent molecule name:", name)
            print ("Current parent molecule path:", path)
            print ("Initial molecule count:", self.__mol_count)

        if name in self.__name_list:
            print ("Molecule %s already found in database at position %d." %(name, self.__name_list.index(name)))
            print ("Not adding molecule to database.")
            return False

        self.__dec_list.append(decompose)
        self.__mol_path.append(path)
        self.__name_list.append(name)

        self.__mol2frag.append(list())
        self.__mol2fragmap.append(list())
        self.__mol2conn.append(list())

        for conn_new_idx, conn_new in enumerate(decompose.connectors):

            found_conn    = False
            found_conn_id = -1

            if self.__conn_count == 0:
                self.__conn2mol.append(list())

                self.__mol2conn[self.__mol_count].append(0)
                self.__conn2mol[0].append(self.__mol_count)

                self.__conn_list.append(conn_new)
                self.__incr_conn()
                continue

            for conn_old_idx, conn_old in  enumerate(self.__conn_list):
                found_conn = are_conn_equal(conn_new, conn_old, self.verbose)
                if found_conn:
                    found_conn_id = conn_old_idx
                    if self.verbose:
                        print ("Current connector", conn_new_idx,)
                        print ("already found at connector id", conn_old_idx)
                    break

            if not found_conn:
                found_conn_id = self.__conn_count
                self.__conn_list.append(conn_new)
                self.__conn2mol.append(list())
                
                self.__incr_conn()

            self.__conn2mol[found_conn_id].append(self.__mol_count)
            self.__mol2conn[self.__mol_count].append(found_conn_id)

        frag_i = 0
        for frag_new, frag_map in zip(decompose.frag_list, decompose.frag_list_map):

            found_frag    = False
            found_frag_id = -1

            ### If we call this routine for the first time
            if self.__frag_count == 0:
                
                self.__frag_list.append(frag_new)

                self.__frag2mol.append(list())
                self.__frag2mol[0].append(self.__mol_count)

                self.__mol2frag[self.__mol_count].append(0)
                self.__mol2fragmap[self.__mol_count].append(frag_map)

                self.__incr_frag()

            ### If we already have stored more than one fragment in the db
            else:
                ### Search fragment
                for frag_id, frag in enumerate(self.__frag_list):
                    found_frag = are_mol_same(frag_new, frag)
                    if found_frag:
                        found_frag_id = frag_id
                        if self.verbose:
                            print ("Current fragment already found at fragment id", frag_id)
                        break

                ### If the fragment is new to the database
                if not found_frag:
                    found_frag_id = self.__frag_count
                    self.__frag_list.append(frag_new)
                    self.__frag2mol.append(list())

                    self.__incr_frag()

                if self.verbose:
                    if not found_frag:
                        print ("Fragment is new to the database and is",)
                        print ("assigned fragment ID", found_frag_id)

                self.__frag2mol[found_frag_id].append(self.__mol_count)
                self.__mol2frag[self.__mol_count].append(found_frag_id)
                self.__mol2fragmap[self.__mol_count].append(frag_map)

            if not found_frag:

                self.__frag2rcapconn.append(list())
                self.__frag2lcapconn.append(list())

                self.__frag2rcapmap.append(list())
                self.__frag2lcapmap.append(list())

                self.__frag2rancmap.append(list())
                self.__frag2lancmap.append(list())

                for cap_conn_idx, cap_list_map, anc_idx_list in zip(decompose.rcap_conn_idx[frag_i],
                                                                    decompose.rcap_list_map[frag_i],
                                                                    decompose.ranc_list[frag_i]):
                    conn_idx = self.__mol2conn[self.__mol_count][cap_conn_idx]
                    self.__frag2rcapconn[found_frag_id].append(conn_idx)
                    self.__frag2rcapmap[found_frag_id].append(cap_list_map)
                    self.__frag2rancmap[found_frag_id].append(anc_idx_list)

                for cap_conn_idx, cap_list_map, anc_idx_list in zip(decompose.lcap_conn_idx[frag_i],
                                                                    decompose.lcap_list_map[frag_i],
                                                                    decompose.lanc_list[frag_i]):
                    conn_idx = self.__mol2conn[self.__mol_count][cap_conn_idx]
                    self.__frag2lcapconn[found_frag_id].append(conn_idx)
                    self.__frag2lcapmap[found_frag_id].append(cap_list_map)
                    self.__frag2lancmap[found_frag_id].append(anc_idx_list)

            frag_i += 1

        for frag2frag_frgs, frag2frag_atms in zip(decompose.frag2frag_frgs,
                                                    decompose.frag2frag_atms):
            frg1 = frag2frag_frgs[0]
            frg2 = frag2frag_frgs[1]
            frag1_idx = self.__mol2frag[-1][frg1]
            frag2_idx = self.__mol2frag[-1][frg2]

            frag1 = self.__frag_list[frag1_idx]
            frag2 = self.__frag_list[frag2_idx]

            frag1_mol = decompose.frag_list[frg1]
            frag2_mol = decompose.frag_list[frg2]

            frag1_matches = frag1.GetSubstructMatches(frag1_mol)[0]
            frag2_matches = frag2.GetSubstructMatches(frag2_mol)[0]

            atms1 = frag1_matches[frag2frag_atms[0]]
            atms2 = frag2_matches[frag2frag_atms[1]]

            if frag1_idx==frag2_idx:
                continue

            self.__frag2frag_frgs.append([frag1_idx, frag2_idx])
            self.__frag2frag_atms.append([atms1, atms2])

        self.__incr_mol()
        return True


class connector(object):

    def __init__(self, name):

        self.__name = name

        self.__conn = None
        self.__surr = None

        self.__lcap       = None
        self.__rcap       = None
        self.__lanc       = list()
        self.__ranc       = list()
        self.__s_lanc     = list()
        self.__r_lanc     = list()
        self.__ring       = True
        self.__terminal   = False
        self.__f_rcap     = list()
        self.__f_lcap     = list()
        self.__terminal_anc = list()

        self.__surr2sur_cap_map = list()

        self.__conn_check      = False
        self.__lcap_check      = False
        self.__rcap_check      = False
        self.__lanc_check      = False
        self.__ranc_check      = False
        self.__ring_check      = False
        self.__terminal_check  = False

        self.__surr_check      = False
        self.__s_lanc_check    = False
        self.__s_ranc_check    = False
        self.__surr_cap_check  = False
        self.__rancmap_check   = False
        self.__lancmap_check   = False
        self.__surr_l_check    = False
        self.__surr_r_check    = False

        self.__surr_cap        = None

        self.__rancmap         = list()
        self.__lancmap         = list()
        self.__rcap_list_check = False
        self.__lcap_list_check = False

    def init_rancmap(self, rancmap):
        if type(rancmap) == list:
            self.__rancmap = rancmap
        else:
            raise TypeError("rancmap must be type list.")

    def init_lancmap(self, lancmap):
        if type(lancmap) == list:
            self.__lancmap = lancmap
        else:
            raise TypeError("lancmap must be type list.")

    def init_surrogate_cap(self):
        
        success = self.__make_surrogate_cap()

        self.__surr_cap_check  = success
        self.__rcap_list_check = success
        self.__lcap_list_check = success
        self.__rancmap_check   = success
        self.__lancmap_check   = success

    def init_surrogate(self, surrogate):
        if type(surrogate) == str:
            self.__surr       = Chem.MolFromSmiles(surrogate)
            self.__surr_check = True
        else:
            raise TypeError("Surrogate must be of type string.")

    def init_s_lanc(self, s_lanc):
        if type(s_lanc) == list:
            self.__s_lanc       = s_lanc
            self.__s_lanc_check = True
        else:
            raise TypeError("s_lanc must be of type string.")

    def init_s_ranc(self, s_ranc):
        if type(s_ranc) == list:
            self.__s_ranc       = s_ranc
            self.__s_ranc_check = True
        else:
            raise TypeError("s_ranc must be of type string.")

    def init_connector(self, connector):
        if type(connector) == str:
            self.__conn       = Chem.MolFromSmiles(connector)
            self.__conn_check = True
        else:
            raise TypeError("Connector must be of type string.")

    def init_lcap(self, lcap):
        if type(lcap) == str:
            self.__lcap       = Chem.MolFromSmiles(lcap)
            self.__lcap_check = True
        else:
            raise TypeError("Lcap must be of type string.")

    def init_rcap(self, rcap):
        if type(rcap) == str:
            self.__rcap       = Chem.MolFromSmiles(rcap)
            self.__rcap_check = True
        else:
            raise TypeError("Rcap must be of type string.")

    def init_lanc(self, lanc):
        if type(lanc) == list:
            self.__lanc       = lanc
            self.__lanc_check = True
        else:
            raise TypeError("Lanc must be of type list.")

    def init_ranc(self, ranc):
        if type(ranc) == list:
            self.__ranc       = ranc
            self.__ranc_check = True
        else:
            raise TypeError("Ranc must be of type list.")

    def init_ring(self, ring):
        if type(ring) == bool:
            self.__ring       = ring
            self.__ring_check = True
        else:
            raise TypeError("ring must be of type bool.")

    def init_terminal(self, terminal):
        if type(terminal) == bool:
            self.__terminal       = terminal
            self.__terminal_check = True
        else:
            raise TypeError("terminal must be of type bool.")

    def init_terminal_anc(self, terminal_anc):
        if type(terminal_anc) == list:
            self.__terminal_anc   = terminal_anc
            self.__terminal_check = True

            self.__rcap_check     = True
            self.__lcap_check     = True

            self.__rcap_list_check = True
            self.__lcap_list_check = True

            self.__s_ranc_check    = True
            self.__s_lanc_check    = True

            self.__surr_cap_check = True
            
            self.__surr_check = True

            self.__surr_r_check = True
            self.__surr_l_check = True

            self.__rcap_list_check = True
            self.__lcap_list_check = True

            self.__rancmap_check = True
            self.__lancmap_check = True
        else:
            raise TypeError("terminal_anc must be of type list.")

    def get_name(self):
        return self.__name

    def get_connector(self):
        return self.__conn

    def get_lcap(self):
        return self.__lcap

    def get_rcap(self):
        return self.__rcap

    def get_lanc(self):
        return self.__lanc

    def get_ranc(self):
        return self.__ranc

    def get_s_lanc(self):
        return self.__s_lanc

    def get_s_ranc(self):
        return self.__s_ranc

    def get_surrogate(self):
        return self.__surr

    def get_ring(self):
        return self.__ring

    def get_terminal(self):
        return self.__terminal

    def get_terminal_anc(self):
        return self.__terminal_anc

    def get_surrogate_cap(self):
        return self.__surr_cap

    def get_rcap_map(self):
        return self.__rcap_list_map

    def get_lcap_map(self):
        return self.__lcap_list_map

    def get_surr2sur_cap_map(self):
        return self.__surr2sur_cap_map

    def get_ranc_map(self):
        return self.__rancmap

    def get_lanc_map(self):
        return self.__lancmap

    def __make_surrogate_cap(self):

        ranc_list = self.get_s_ranc()
        lanc_list = self.get_s_lanc()
        rcap_mol  = self.get_rcap()
        lcap_mol  = self.get_lcap()
        conn_mol  = self.get_surrogate()

        ranc_conn = list()
        lanc_conn = list()

        ranc_bond = list()
        lanc_bond = list()

        ranc_bondtype = list()
        lanc_bondtype = list()

        ranc_parity = list()
        lanc_parity = list()

        ranc_chirality = list()
        lanc_chirality = list()

        for atm_idx in ranc_list:
            atm = conn_mol.GetAtomWithIdx(atm_idx)
            for neighbor in atm.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in ranc_list:
                    bond = conn_mol.GetBondBetweenAtoms(atm_idx,\
                                                        neighbor_idx)

                    ranc_conn.append(neighbor_idx)
                    ranc_bondtype.append(bond.GetBondType())
                    ranc_parity.append(get_parity(conn_mol, neighbor_idx))
                    ranc_chirality.append(neighbor.GetChiralTag())
                    ranc_bond.append(bond.GetIdx())

        for atm_idx in lanc_list:
            atm = conn_mol.GetAtomWithIdx(atm_idx)
            for neighbor in atm.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in lanc_list:
                    bond = conn_mol.GetBondBetweenAtoms(atm_idx,\
                                                        neighbor_idx)

                    lanc_conn.append(neighbor_idx)
                    lanc_bondtype.append(bond.GetBondType())
                    lanc_parity.append(get_parity(conn_mol, neighbor_idx))
                    lanc_chirality.append(neighbor.GetChiralTag())
                    lanc_bond.append(bond.GetIdx())

        frags = Chem.FragmentOnBonds(conn_mol,
                                     ranc_bond+lanc_bond,
                                     addDummies=False)

        frags_mol, mol_frags_atm_idxs = Chem.GetMolFrags(frags, asMols=True),\
                                        Chem.GetMolFrags(frags, asMols=False)

        rcap_length = rcap_mol.GetNumAtoms()
        lcap_length = lcap_mol.GetNumAtoms()

        ###
        ### Loop over all fragments generated by FragmentOnBonds
        ###
        ### frag_mol    : rdkit Mol instance
        ### mol_atm_idxs: atom indices of the original molecules that
        ###               generated the molecule
        ###

        for frag_mol, mol_atm_idxs in zip(frags_mol, mol_frags_atm_idxs):

            frag_length = frag_mol.GetNumAtoms()
            rcap_start  = frag_length
            lcap_start  = frag_length + rcap_length

            combo_mol   = Chem.CombineMols(frag_mol, rcap_mol)
            combo_mol   = Chem.CombineMols(combo_mol, lcap_mol)
            e_combo_mol = Chem.RWMol(combo_mol)

            is_frag    = False

            ranc_check = list()
            lanc_check = list()

            rcap_list_map = range(frag_length,
                                  frag_length+rcap_length)
            lcap_list_map = range(frag_length+rcap_length, 
                                  frag_length+rcap_length+lcap_length)

            ###
            ### Loop over all atoms in the fragment
            ###
            ### atm_idx    : atom index in the fragment indexing scheme
            ### mol_atm_idx: atom index in the original molecule indexing scheme
            ###

            for atm_idx, mol_atm_idx in enumerate(mol_atm_idxs):

                if mol_atm_idx in ranc_conn:
                    idx = ranc_conn.index(mol_atm_idx)
                    e_combo_mol.AddBond(atm_idx,
                                        rcap_start,
                                        ranc_bondtype[idx])

                    ranc_check.append([atm_idx, mol_atm_idx, ranc_parity[idx]])

                    is_frag = True
                if mol_atm_idx in lanc_conn:
                    idx = lanc_conn.index(mol_atm_idx)
                    e_combo_mol.AddBond(atm_idx,
                                        lcap_start,
                                        lanc_bondtype[idx])

                    lanc_check.append([atm_idx, mol_atm_idx, lanc_parity[idx]])
                    is_frag = True

            combo_mol = e_combo_mol.GetMol()

            if is_frag:

                for atm_idx, mol_atm_idx, parity in ranc_check:
                    if get_parity(combo_mol, atm_idx) != parity:
                        atm        = combo_mol.GetAtomWithIdx(atm_idx)
                        new_chiral = None
                        if ranc_chirality == CHI_TETRAHEDRAL_CCW:
                            new_chiral = CHI_TETRAHEDRAL_CW
                        elif ranc_chirality == CHI_TETRAHEDRAL_CW:
                            new_chiral = CHI_TETRAHEDRAL_CCW
                        if new_chiral != None:
                            atm.SetChiralTag(new_chiral)

                for atm_idx, mol_atm_idx, parity in lanc_check:
                    if get_parity(combo_mol, atm_idx) != parity:
                        atm        = combo_mol.GetAtomWithIdx(atm_idx)
                        new_chiral = None
                        if lanc_chirality == CHI_TETRAHEDRAL_CCW:
                            new_chiral = CHI_TETRAHEDRAL_CW
                        elif lanc_chirality == CHI_TETRAHEDRAL_CW:
                            new_chiral = CHI_TETRAHEDRAL_CCW
                        if new_chiral != None:
                            atm.SetChiralTag(new_chiral)

                if get_parity(combo_mol, rcap_start) != get_parity(rcap_mol, 0):
                    atm        = combo_mol.GetAtomWithIdx(rcap_start)
                    chirality  = atm.GetChiralTag()
                    new_chiral = None
                    if chirality == CHI_TETRAHEDRAL_CCW:
                        new_chiral = CHI_TETRAHEDRAL_CW
                    elif chirality == CHI_TETRAHEDRAL_CW:
                        new_chiral = CHI_TETRAHEDRAL_CCW
                    if new_chiral != None:
                        atm.SetChiralTag(new_chiral)

                if get_parity(combo_mol, lcap_start) != get_parity(lcap_mol, 0):
                    atm        = combo_mol.GetAtomWithIdx(lcap_start)
                    chirality  = atm.GetChiralTag()
                    new_chiral = None
                    if chirality == CHI_TETRAHEDRAL_CCW:
                        new_chiral = CHI_TETRAHEDRAL_CW
                    elif chirality == CHI_TETRAHEDRAL_CW:
                        new_chiral = CHI_TETRAHEDRAL_CCW
                    if new_chiral != None:
                        atm.SetChiralTag(new_chiral)

                self.__surr_cap = copy.copy(combo_mol)

                self.__rcap_list_map = rcap_list_map
                self.__lcap_list_map = lcap_list_map

                self.__surr2sur_cap_map = mol_atm_idxs

                return True

        return False


    def check(self, verbose=False):

        anc_align = True

        if self.__conn_check \
        and self.__ranc_check \
        and self.__lanc_check:
            natms       = self.__conn.GetNumAtoms()
            n_ranc_atms = len(self.__ranc )
            n_lanc_atms = len(self.__lanc)
            if natms != (n_ranc_atms+n_lanc_atms):
                anc_align = False
            elif max(self.__ranc) > natms: anc_align = False
            elif max(self.__lanc) > natms: anc_align = False
            else:
                for i in range(natms):
                    if i in self.__ranc \
                    and i in self.__lanc:
                        anc_align = False
                    if not (i in self.__ranc \
                    or i in self.__lanc):
                        anc_align = False

        if verbose:
            print ("__conn_check:", self.__conn_check)
            print ("__lcap_check:", self.__lcap_check)
            print ("__rcap_check:", self.__rcap_check)
            print ("__lanc_check:", self.__lanc_check)
            print ("__ranc_check:", self.__ranc_check)
            print ("__ring_check:", self.__ring_check)
            print ("__terminal_check:", self.__terminal_check)
            print ("__s_lanc_check:", self.__s_lanc_check)
            print ("__s_ranc_check:", self.__s_ranc_check)
            print ("__surr_check:", self.__surr_check)
            print ("__surr_cap_check:", self.__surr_cap_check)
            print ("__rcap_list_check:", self.__rcap_list_check)
            print ("__lcap_list_check:", self.__lcap_list_check)
            print ("__lancmap_check:", self.__lancmap_check)
            print ("__rancmap_check:", self.__rancmap_check)
            print ("anc_align:", anc_align)

        if self.__conn_check \
           and self.__lcap_check \
           and self.__rcap_check \
           and self.__lanc_check \
           and self.__ranc_check \
           and self.__ring_check \
           and self.__terminal_check \
           and self.__s_lanc_check \
           and self.__s_ranc_check \
           and self.__surr_check \
           and self.__surr_cap_check \
           and self.__rcap_list_check \
           and self.__lcap_list_check \
           and self.__rancmap_check \
           and self.__lancmap_check \
           and anc_align:
            return True

        else:
            return False