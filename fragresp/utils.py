from rdkit import Chem
import os
from subprocess import call
from collections import OrderedDict

from fragresp.aux_progs import ante_exe
from fragresp.aux_progs import respgen_exe
from fragresp.aux_progs import espgen_exe
from fragresp.aux_progs import resp_exe

from molvs import tautomer
from molvs import standardize
from molvs import standardize_smiles

from rdkit import Chem
from rdkit.Chem import AllChem

class mol2_reader(object):

    def __init__(self, Path):

        self.path    = Path
        
        self.N_atoms = 0
        
        self.atomid   = list()
        self.atomname = list()
        self.crds     = list()
        self.types    = list()
        self.resid    = list()
        self.resname  = list()
        self.charges  = list()

        self.__infile = None

        self.__read()

    def __read(self):

        found_start   = False
        self.__infile = open(self.path, "r")

        for l in self.__infile:
            line = l.rstrip().split()
            if len(line) == 0:
                continue
            if not found_start and line[0] == "@<TRIPOS>ATOM":
                found_start = True
                continue
            if line[0] == "@<TRIPOS>BOND":
                break
            if found_start:
                self.atomid.append(int(line[0]))
                self.atomname.append(str(line[1]))
                self.crds.append([float(line[2]), float(line[3]), float(line[4])])
                self.types.append(str(line[5]))
                self.resid.append(int(line[6]))
                self.resname.append(str(line[7]))
                self.charges.append(float(line[8]))

        self.N_atoms = len(self.atomid)


    def write(self, path):

        found_start = False
        found_end   = False

        self.__infile.seek(0)

        with open(path, "w") as outfile:
            for l in self.__infile:
                line = l.rstrip().split()
                if len(line)>0 and line[0] == "@<TRIPOS>BOND":
                    found_end = True
                if found_start and not found_end:
                    continue
                if len(line)>0 and not found_start and line[0] == "@<TRIPOS>ATOM":
                    outfile.write("@<TRIPOS>ATOM")
                    outfile.write("\n")
                    for i in range(self.N_atoms):
                        outfile.write(" %d " %self.atomid[i])
                        outfile.write(" %s " %self.atomname[i])
                        outfile.write(" %f " %self.crds[i][0])
                        outfile.write(" %f " %self.crds[i][1])
                        outfile.write(" %f " %self.crds[i][2])
                        outfile.write(" %s " %self.types[i])
                        outfile.write(" %d " %self.resid[i])
                        outfile.write(" %s " %self.resname[i])
                        outfile.write(" %f " %self.charges[i])
                        outfile.write("\n")

                    found_start = True
                    continue

                outfile.write(l)


class logger(object):

    def __init__(self, name):

        self.name     = name
        self.__isopen = False
        self.__open()

    def __open(self):

        self.__file   = open(self.name, "w")
        self.__isopen = True

    def close(self):

        self.__file.close()
        self.__isopen = False

    def log(self, text=''):

        if self.__isopen:

            self.__file.write(text)
            self.__file.write('\n')

        else:

            raise IOError("logger object is already closed.")


def get_frag_bonds(Mol, Anc_list, Connector_list, Ring, Verbose=False):

    ### This function takes a list of ancor and connector atoms
    ### and returns 
    ### A) a list of the bond indices of the bonds to be broken
    ### B) a list of the atom indices of the atoms in the fragment
    ###    that are part of the broken bond.
    ### Note, that it can process only one ancor/connector pair at
    ### a time.

    ### Mol            : instance of rdkit molecule object
    ### Anc_list       : list of ancor atom indices
    ### Connector_list : list of connector atom indices
    ### Ring           : boolean indicating if in-ring bonds must be preserved
    ###
    ### returns tuple ( bond_list, atom_list )
    ###
    ### bond_list: bond indices of bonds between anchors and connectors
    ### atom_list: tuple of atom indices of anchor atom and connector atom:
    ###            (anchor_atom, connector_atom)

    bond_list = list()
    atom_list = list()

    for anc_idx in Anc_list:
        anc           = Mol.GetAtomWithIdx(anc_idx)
        anc_neighbors = anc.GetNeighbors()
        if Verbose:
            print ("List of neighbor atoms at anchor atom", anc_idx,":",)
            for anc_neighbor in anc_neighbors:
                print (anc_neighbor.GetIdx(),)
            print ("")
        ### Loop over all neighboring atoms of anchor atoms
        for anc_neighbor in anc_neighbors:
            anc_neighbor_idx = anc_neighbor.GetIdx()
            if anc_neighbor_idx in Connector_list \
            and anc_neighbor_idx not in Anc_list:
                bond = Mol.GetBondBetweenAtoms(anc_idx, anc_neighbor_idx)
                if Ring:
                    if not bond.IsInRing():
                        bond_list.append(bond.GetIdx())
                        atom_list.append([anc_idx, anc_neighbor_idx])
                    #if not (anc.IsInRing() and anc_neighbor.IsInRing()):
                    #    bond_list.append(Mol.GetBondBetweenAtoms(anc_idx, anc_neighbor_idx).GetIdx())
                    #    atom_list.append([anc_idx, anc_neighbor_idx])
                    elif Verbose:
                        print ("Bond between atoms", anc_idx, "and", anc_neighbor_idx, "not broken",)
                        print ("due to ring bond preservation rule.")
                else:
                    #bond_list.append(Mol.GetBondBetweenAtoms(anc_idx, anc_neighbor_idx).GetIdx())
                    #atom_list.append([anc_idx, anc_neighbor_idx])
                    bond_list.append(bond.GetIdx())
                    atom_list.append([anc_idx, anc_neighbor_idx])

    return bond_list, atom_list


def are_conn_equal(conn1, conn2, verbose=False):

    if not conn1.check(verbose):
        raise Warning ("Cannot compare connectors. Conn1 is not fully defined.")
        return False
    if not conn2.check(verbose):
        raise Warning ("Cannot compare connectors. Conn2 is not fully defined.")
        return False

    if not are_mol_same(conn1.get_connector(), conn2.get_connector()):
        if verbose:
            "Connectors not same because different connector motifs."
        return False

    if not conn1.get_terminal() and not conn1.get_terminal():
        if not are_mol_same(conn1.get_lcap(), conn2.get_lcap()):
            if verbose:
                "Connectors not same because different Lcaps."
            return False

        if not are_mol_same(conn1.get_rcap(), conn2.get_rcap()):
            if verbose:
                "Connectors not same because different Rcaps."
            return False

        if not are_mol_same(conn1.get_surrogate(), conn2.get_surrogate()):
            if verbose:
                "Connectors not same because different surrogates."
            return False

        if not are_mol_same(conn1.get_surrogate_cap(), conn2.get_surrogate_cap()):
            if verbose:
                "Connectors not same because different capped surrogates."
            return False

    for c1 in conn1.get_lanc():
        if c1 not in conn2.get_lanc():
            if verbose:
                "Connectors not same because different L anchors."
            return False
    for c2 in conn2.get_lanc():
        if c2 not in conn1.get_lanc():
            if verbose:
                "Connectors not same because different L anchors."
            return False

    for c1 in conn1.get_ranc():
        if c1 not in conn2.get_ranc():
            if verbose:
                "Connectors not same because different R anchors."
            return False
    for c2 in conn2.get_ranc():
        if c2 not in conn1.get_ranc():
            if verbose:
                "Connectors not same because different R anchors."
            return False

    if conn1.get_ring() != conn2.get_ring():
        if verbose:
            "Connectors not same because different ring preservation rule."
        return False

    if conn1.get_terminal() != conn2.get_terminal():
        if verbose:
            "Connectors not same because different treatment of terminal linkers."
        return False

    for c1 in conn1.get_terminal_anc():
        if c1 not in conn2.get_terminal_anc():
            if verbose:
                "Connectors not same because different L terminal anchors."
            return False

    return True


def are_mol_same(mol1, mol2):

    is_same     = False
    mol1_length = mol1.GetNumAtoms()
    mol2_length = mol2.GetNumAtoms()
    mol1_atoms  = range(mol1_length)
    ### If both have same number of atoms.
    if mol1_length == mol2_length:
        matches = mol1.GetSubstructMatches(mol2, useChirality=True)
        ### If there is only one substructure match, which also has
        ### as many atoms as the other molecule.
        if len(matches) == 1 and len(matches[0]) == mol1_length:
            is_same = True
    return is_same


def get_parity(Mol, Atm_idx):

    canonical_rank = list()
    neighbor_list  = list()
    neighbor_rank  = list()
    string_rank    = list(Chem.CanonicalRankAtoms(Mol, breakTies=False))
    for rank in string_rank:
        canonical_rank.append(int(rank))
    del string_rank

    for bond in Mol.GetAtomWithIdx(Atm_idx).GetBonds():
        neighbor_idx = bond.GetOtherAtomIdx(Atm_idx)
        neighbor_list.append(neighbor_idx)
        neighbor_rank.append(canonical_rank[neighbor_idx])

    ### See also http://www.dalkescientific.com/writings/diary/archive/2016/08/14/fragment_chiral_molecules.html
    N         = len(neighbor_rank)
    num_swaps = 0
    for i in range(N-1):
        for j in range(i+1, N):
            if neighbor_rank[i] > neighbor_rank[j]:
                neighbor_rank[i], neighbor_rank[j] = neighbor_rank[j], neighbor_rank[i]
                num_swaps += 1
    return num_swaps % 2


def cat_esp(esp_list, target_esp):

    with open(target_esp, "w") as main_esp:
        for esp_file in esp_list:
            with open(esp_file) as new_esp:
                for line in new_esp:
                    main_esp.write(line)

    return True


def preprocess_esp(log_path, stdout=None, stderr=None):

    if stdout==None:
        stdout=open(os.devnull, 'w')
    if stderr==None:
        stderr=open(os.devnull, 'w')

    mainpath          = os.path.splitext(log_path)[0]
    ante_args_general = ["-fi", "gout", "-pf", "y", "-at", "sybyl", "-dr", "no", "-j", "5"]

    mol2_path   = mainpath+".mol2"
    esp_path    = mainpath+".esp"
    ante_args   = [ante_exe, "-i", log_path, "-o", mol2_path, "-fo", "mol2"] + ante_args_general
    espgen_args = [espgen_exe, "-i", log_path, "-o", esp_path]

    call(ante_args, stdout=stdout, stderr=stderr)
    call(espgen_args, stdout=stdout, stderr=stderr)

    return esp_path, mol2_path


def fix_groups(mol, verbose=False):

    """
    Return fixed copy of mol

    fix_dict: key=name of the fix, value=[wrong mol, correct mol]
    """

    fix_dict = {
    'nitro' : [Chem.MolFromSmiles("N([O-])[O-]"), Chem.MolFromSmiles("[N+](=O)[O-]")]
    }
    for name, mol_replace in fix_dict.items():
        if mol.HasSubstructMatch(mol_replace[0]):
            if verbose:
                print ("Found bad %s substructure. Attempting correction." %name)
            new = AllChem.ReplaceSubstructs(mol, mol_replace[0], mol_replace[1])
            mol = new[0]
    return mol


def canonicalize_tautomers_old(rank_list, mol):

    for atm_idx in range(mol.GetNumAtoms()):
        atm     = mol.GetAtomWithIdx(atm_idx)
        atm_num = atm.GetAtomicNum()
        if atm_num == 1:
            continue
        nghbr_cnt = 0
        for nghbr in atm.GetNeighbors():
            nghbr_num = nghbr.GetAtomicNum()
            if nghbr_num > 1:
                nghbr_cnt += 1
        if nghbr_cnt < 2:
            continue
        atm_chg = atm.GetFormalCharge()
        atm_hyb = atm.GetHybridization()

        idx = list()
        num = list()
        chg = list()
        hyb = list()
        cnt = 0
        for nghbr in atm.GetNeighbors():
            idx.append(nghbr.GetIdx())
            num.append(nghbr.GetAtomicNum())
            chg.append(nghbr.GetFormalCharge())
            hyb.append(nghbr.GetHybridization())
            cnt += 1

        ### Find neighbor atoms 'nghbr' that have same hybridization as 
        ### parent atom 'atm'
        for i in range(cnt):
            if hyb[i] == atm_hyb:
                for j in range(cnt):
                    if j > i \
                    and hyb[j] == hyb[i] \
                    and num[j] == num[i]:
                        rank_list[idx[j]] = rank_list[idx[i]]


def canonicalize_tautomers(rank_list, mol):

    canon = tautomer.TautomerCanonicalizer()
    mol_t = canon.canonicalize(mol)
    rank_list = list()
    string_rank = list(Chem.CanonicalRankAtoms(mol_t, breakTies=False))
    for rank in string_rank:
        rank_list.append(int(rank))
    del string_rank