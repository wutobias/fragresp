import os
from subprocess import call
import fortranformat
import numpy as np
from rdkit import Chem

from fragresp.aux_progs import ante_exe
from fragresp.aux_progs import respgen_exe
from fragresp.aux_progs import espgen_exe
from fragresp.aux_progs import resp_exe
from fragresp.utils import cat_esp
from fragresp.utils import preprocess_esp
from fragresp.utils import canonicalize_tautomers
from fragresp.utils import fix_groups

class write_resp_in(object):

    def __init__(self):

        ### Title of the resp input file. *Not* file name.
        self._title             = "Resp input file generated with FragResp"

        ### list of length N (=number of molecules)
        self._mol_list          = list()
        self._mol_weight_list   = list()
        self._mol_name_list     = list()
        self._mol_charge_list   = list()

        ### list of constrain groups
        ### _group_atom_list contains atom indices of atoms that are
        ### bound to groups. The list _group_mol_list has same length
        ### as _group_atom_list and stores the corresponding mol_id
        ### of every atom in _group_atom_list. The list _group_id_list
        ### contains the uniq group id for every atom in _group_atom_list.
        ### The list _group_charge_list contains the total charge 
        ### of every group. This list has length self._group_count.
        self._group_atom_list   = list()
        self._group_mol_list    = list()
        self._group_charge_list = list()
        self._group_id_list     = list()

        ### The list _atom_equiv_list_i contains indices of atoms in a group i
        ### that should be equal to other atoms in a group j, which are stored
        ### correspondingly in _atom_equiv_mol_j. These atom indices are assigned
        ### to molecules, whose indices are stored in _atom_equiv_mol_i and 
        ### _atom_equiv_mol_j, respectively.
        self._atom_equiv_list_i = list()
        self._atom_equiv_list_j = list()
        self._mol_equiv_mol_i   = list()
        self._mol_equiv_mol_j   = list()

        self._mol_count         = 0
        self._group_count       = 0
        self._atom_equiv_count  = 0

        self._intermol1         = list()
        self._intermol2         = list()

        self.stage              = 'resp1'

        self.groups_frozen      = False
        self.noh_frozen         = False
        self.h_equiv            = False
        self.all_equiv          = False
        self.unfreeze_all       = False

        ### The self._free_list contains a list (can also be empty list) for each
        ### molecule in the list self._free_list_mol, holding the indices of the 
        ### atoms that are not allowed to be freezed or constraint in any way.
        self._free_list     = list()
        self._free_list_mol = list()

    def add_free(self, mol_i, free):

        if type(mol_i) != int:
            raise TypeError("mol_i must be of type int")
        if type(free) != list:
            raise TypeError("free must be of type list")

        if mol_i in self._free_list_mol:
            mol_i_idx = self._free_list_mol.index(mol_i)
            self._free_list[mol_i_idx] += free
        else:
            self._free_list_mol.append(mol_i)
            self._free_list.append(free)

    def add_mol(self, mol, charge, weight=1.0, name=None):
        self._mol_list.append(mol)
        self._mol_weight_list.append(weight)
        self._mol_charge_list.append(charge)
        if name == None:
            name = "Organic Molecule %d" %self._mol_count
        self._mol_name_list.append(name)
        self._mol_count += 1

    def add_group(self, atom_i_list, mol_i_list, charge):
        """
        atom_i_list: List of atom indices that should
                     go into one group.
        mol_i_list : List of mol indices of atoms in atom_i_list
                     that should go into one group
        charge     : Charge of that group
        """

        if type(atom_i_list) != list:
            raise TypeError("atom_i_list must be of type list but is of type %s" %type(atom_i_list))

        if type(mol_i_list) != list:
            raise TypeError("mol_i_list must be of type list but is of type %s" %type(mol_i_list))

        if len(atom_i_list) != len(mol_i_list):
            raise Exception("atom_i_list (len %d) and mol_i_list (len %d) must be of same legnth." \
                          %(len(atom_i_list), len(mol_i_list)))

        if type(charge) != float:
            raise TypeError("charge must be of type float but is of type %s" %type(charge))

        self._group_atom_list.extend(atom_i_list)
        self._group_mol_list.extend(mol_i_list)
        self._group_id_list.extend([self._group_count]*len(atom_i_list))
        self._group_charge_list.append(charge)
        self._group_count += 1

    def add_intermolecular(self, mol_1, mol_2):

        if mol_1 in self._intermol1:
            index = self._intermol1.index(mol_1)
            if mol_2 not in self._intermol2[index]:
                self._intermol2[index].append(mol_2)
        else:
            self._intermol1.append(mol_1)
            self._intermol2.append(list())
            self._intermol2[-1].append(mol_2)

    def add_atom_equiv(self, atom_list_i, atom_list_j, mol_list_i, mol_list_j):
        
        if type(atom_list_i) != list:
            raise TypeError("atom_list_i must be of type list but is of type %s." %type(atom_list_i))

        if type(atom_list_j) != list:
            raise TypeError("atom_list_j must be of type list but is of type %s." %type(atom_list_j))

        if type(mol_list_i) != list:
            raise TypeError("mol_list_i must be of type list but is of type %s." %type(mol_list_i))

        if type(mol_list_j) != list:
            raise TypeError("mol_list_j must be of type list but is of type %s." %type(mol_list_j))

        if len(atom_list_i) != len(atom_list_j):
            raise ValueError("atom_list_i is of length %d and atom_list_j is of length %d. They should be equal."\
                              %(len(atom_list_i), len(atom_list_j)))

        if len(mol_list_i) != len(mol_list_j):
            raise ValueError("mol_list_i is of length %d and mol_list_j is of length %d. They should be equal."\
                              %(len(mol_list_i), len(mol_list_j)))

        if len(atom_list_i) != len(mol_list_i):
            raise ValueError("atom_list_i is of length %d and mol_list_i is of length %d. They should be equal."\
                              %(len(atom_list_i), len(mol_list_i)))

        if len(atom_list_j) != len(mol_list_j):
            raise ValueError("atom_list_j is of length %d and mol_list_j is of length %d. They should be equal."\
                              %(len(atom_list_j), len(mol_list_j)))

        self._atom_equiv_list_i.extend(atom_list_i)
        self._atom_equiv_list_j.extend(atom_list_j)
        self._mol_equiv_mol_i.extend(mol_list_i)
        self._mol_equiv_mol_j.extend(mol_list_j)

        self._atom_equiv_count += len(atom_list_i)

    def set_title(self, title_string):
        self._title = title_string

    def set_stage(self, stage):

        valid_choices = ['resp1_surr',
                         'resp2_surr',
                         'resp1',
                         'resp2']

        if type(stage) != str:
            raise ValueError("stage must be of type string but is of type %s" %type(stage))

        if stage not in valid_choices:
            raise IOError("stage is %s, but must be one of %s" %(stage, valid_choices))

        self.stage = stage

    def get_title(self):

        """
        Return title.
        """

        return self._title

    def get_mol(self):

        """
        Return section containing element types, fitting weight, molecule
        title, number of atoms and atom equivalencing.

        groups_frozen: Freeze charges in groups to the values in qin file,
                       typcially obtained from previous resp run.

        h_equiv      : Fit charges of degenerate hydrogen atoms together

        all_equiv    : Freeze charges of all degenerate atoms together. If
                       this is activated, and h_equiv is deactivated, only
                       heavy-atom atomic centers will be fitted together.
        """

        line_2I5   = fortranformat.FortranRecordWriter('2I5')

        _tmp_str   = list()

        for mol_i in range(self._mol_count):
            mol = self._mol_list[mol_i]

            _tmp_str.append('  %f\n' %self._mol_weight_list[mol_i])
            _tmp_str.append('  %s\n' %self._mol_name_list[mol_i])

            _charge    = float(self._mol_charge_list[mol_i])
            _charge    = round(_charge)
            _charge    = int(_charge)
            _natoms    = mol.GetNumAtoms()
            _tmp_str.append(line_2I5.write([_charge, _natoms]))
            _tmp_str.append('\n')

            canonical_rank = list()
            string_rank    = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
            for rank in string_rank:
                canonical_rank.append(int(rank))
            del string_rank
            ### This really never worked perfectly...
            canonicalize_tautomers(canonical_rank, mol)

            index_list = np.arange(_natoms)

            if mol_i not in self._intermol1:
                for atom_i in index_list:
                    atom   = mol.GetAtomWithIdx(int(atom_i))
                    at_num = atom.GetAtomicNum()
                    _tmp_str.append(line_2I5.write([at_num, 0]))
                    _tmp_str.append('\n')

            else:
                for atom_i in index_list:
                    atom   = mol.GetAtomWithIdx(int(atom_i))
                    at_num = atom.GetAtomicNum()
                    if mol_i in self._free_list_mol:
                        mol_i_idx = self._free_list_mol.index(mol_i)
                        if atom_i in self._free_list[mol_i_idx]:
                            _tmp_str.append(line_2I5.write([at_num, 0]))
                            _tmp_str.append('\n')
                            continue

                    placed_frozen = False

                    if self.unfreeze_all:
                        _tmp_str.append(line_2I5.write([at_num, 0]))

                    elif self.groups_frozen:
                        ### Check if atom itself is in group
                        for index, atom_j in enumerate(self._group_atom_list):
                            if atom_j==atom_i \
                            and self._group_mol_list[index] == mol_i:
                                if self.noh_frozen and at_num != 1:
                                    _tmp_str.append(line_2I5.write([at_num, -1]))
                                elif not self.h_groups_frozen and at_num == 1:
                                    if self.h_equiv:
                                        canon_eq_bool = np.isin(canonical_rank, canonical_rank[atom_i])
                                        canon_eq_int  = index_list[canon_eq_bool]
                                        if atom_i == canon_eq_int[0]:
                                            _tmp_str.append(line_2I5.write([at_num, 0]))
                                        else:
                                            _tmp_str.append(line_2I5.write([at_num, canon_eq_int[0]+1]))
                                    else:
                                        _tmp_str.append(line_2I5.write([at_num, 0]))
                                else:
                                    _tmp_str.append(line_2I5.write([at_num, 0]))
                                placed_frozen = True
                                break

                    if not placed_frozen:
                        if self.noh_frozen and at_num != 1:
                            _tmp_str.append(line_2I5.write([at_num, -1]))
                        elif (self.h_equiv and at_num == 1) \
                        or (self.all_equiv and not self.h_equiv and at_num != 1):
                            ### canon_eq_bool is True for all atoms that are canonically
                            ### equal to atom_i (including atom_i itself).
                            ### canon_eq_int holds atom indices of all atoms that are
                            ### canonically equal to atom_i (including atom_i itself).
                            canon_eq_bool = np.isin(canonical_rank, canonical_rank[atom_i])
                            canon_eq_int  = index_list[canon_eq_bool]
                            ### This is fulfilled only when we encounter this canoncial
                            ### rank (stored in canonical_rank[atom_i]) for the first
                            ### time in this molecule. It will tell resp to let that 
                            ### atom center vary independly.
                            if atom_i == canon_eq_int[0]:
                                _tmp_str.append(line_2I5.write([at_num, 0]))
                            ### If current atom atom_i is equivalent to another atom
                            ### which is present in a different group than atom_i, then
                            ### we should not equivalence constraints on these two atoms.
                            ### If we already have encountered this canoncial rank before
                            ### freeze atom_i to the atom that was our first encounter with
                            ### this canonical rank. Note, that resp expects atom counting 
                            ### to start at 1, *not* 0.
                            else:
                                _tmp_str.append(line_2I5.write([at_num, canon_eq_int[0]+1]))
#                            else:
#                                found_in_group = False
#                                for index2, atom_j in enumerate(self._group_atom_list):
#                                    if canon_eq_int[0]==atom_j \
#                                    and self._group_mol_list[index2]==mol_i:
#                                        for index, atom_k in enumerate(self._group_atom_list):
#                                            if atom_i==atom_k \
#                                            and self._group_mol_list[index]==mol_i:
#                                                if index2 == index:
#                                                    _tmp_str.append(line_2I5.write([at_num, canon_eq_int[0]+1]))
#                                                else:
#                                                    _tmp_str.append(line_2I5.write([at_num, 0]))
#                                                found_in_group = True
#                                            if found_in_group:
#                                                break
#                                    if found_in_group:
#                                        break
#
#                                if not found_in_group:
#                                    _tmp_str.append(line_2I5.write([at_num, 0]))

                        else:
                            _tmp_str.append(line_2I5.write([at_num, 0]))

                    _tmp_str.append('\n')

            if self._mol_count>1:
                _tmp_str.append('\n')

        return ''.join(_tmp_str)

    def get_group(self):

        line_I5_F105     = fortranformat.FortranRecordWriter('I5,F10.5')
        line_16I5        = fortranformat.FortranRecordWriter('16I5')

        _tmp_str = list()

        group_index_list = np.arange(len(self._group_atom_list))

        for group_i in range(self._group_count):

            ### group_valids_bool is True whereever an atom is in 
            ### the group with group id group_i.
            ### group_valids_int gives the list indices of every atom
            ### in self._group_atom_list that is in group group_i.
            group_valids_bool = np.isin(self._group_id_list, group_i)
            group_valids_int  = group_index_list[group_valids_bool]

            group_size   = group_valids_int.shape[0]
            group_charge = self._group_charge_list[group_i]

            _tmp_str.append(line_I5_F105.write([group_size, group_charge]))
            _tmp_str.append('\n')

            atom_list = list()

            for valids_i in group_valids_int:

                atom_i = self._group_atom_list[valids_i]
                mol_id = self._group_mol_list[valids_i]

                atom_list.append(mol_id+1)
                atom_list.append(atom_i+1)

            _tmp_str.append(line_16I5.write(atom_list))

            _tmp_str.append('\n')

        _tmp_str.append('\n')

        return ''.join(_tmp_str)

    def get_intermolecular(self):

        line_16I5 = fortranformat.FortranRecordWriter('16I5')
        line_I5   = fortranformat.FortranRecordWriter('I5')

        _tmp_str  = list()

        for mol_1, mol_2_list in zip(self._intermol1,\
                                     self._intermol2):
            ### First, check if the molecule has more
            ### than one conformer. If not, continue with
            ### for-loop.
            N_mols = 1
            for mol_2 in mol_2_list:
                if mol_2 != mol_1:
                    N_mols += 1
            if N_mols == 1:
                continue
            for atm_idx in range(self._mol_list[mol_1].GetNumAtoms()):
                _mols_list = list()
                _mols_list.append(mol_1+1)
                _mols_list.append(atm_idx+1)
                for mol_2 in mol_2_list:
                    if mol_2 != mol_1:
                        _mols_list.append(mol_2+1)
                        _mols_list.append(atm_idx+1)

                if len(_mols_list)> 0:
                    _tmp_str.append(line_I5.write([N_mols]))
                    _tmp_str.append('\n')
                    _tmp_str.append(line_16I5.write(_mols_list))
                _tmp_str.append('\n')

        return ''.join(_tmp_str)

    def get_atom_equiv(self):

        line_16I5 = fortranformat.FortranRecordWriter('16I5')
        line_I5   = fortranformat.FortranRecordWriter('I5')

        _tmp_str  = list()

        mol_atom_i = [x for x in zip(self._mol_equiv_mol_i, self._atom_equiv_list_i)]
        mol_atom_j = [x for x in zip(self._mol_equiv_mol_j, self._atom_equiv_list_j)]

        exclude    = list()

        for i in range(self._atom_equiv_count):

            if i in exclude:
                continue

            N_mols = 0

            _mols_list = list()
            _mols_list.append(mol_atom_i[i][0]+1)## 9   1
            _mols_list.append(mol_atom_i[i][1]+1)
            _mols_list.append(mol_atom_j[i][0]+1)## 4   7
            _mols_list.append(mol_atom_j[i][1]+1)

            N_mols += 2

            exclude.append(i)

            for j in range(self._atom_equiv_count):

                if j in exclude:
                    continue

                if mol_atom_j[j] == mol_atom_i[i] \
                and mol_atom_i[j] != mol_atom_i[i]:
                    _mols_list.append(mol_atom_i[j][0]+1)
                    _mols_list.append(mol_atom_i[j][1]+1)

                    N_mols += 1

                    exclude.append(j)

                if mol_atom_i[j] == mol_atom_i[i] \
                and mol_atom_j[j] != mol_atom_i[i]:
                    _mols_list.append(mol_atom_j[j][0]+1)
                    _mols_list.append(mol_atom_j[j][1]+1)

                    N_mols += 1

                    exclude.append(j)

            if len(_mols_list) > 0:
                _tmp_str.append(line_I5.write([N_mols]))
                _tmp_str.append('\n')
                _tmp_str.append(line_16I5.write(_mols_list))
            _tmp_str.append('\n')

        return ''.join(_tmp_str)

    def get_cntrl(self):

        _tmp_str = list()
        _tmp_str.append(' &cntrl')
        _tmp_str.append('\n')
        _tmp_str.append('\n')

        ### self.groups_frozen  : Atoms in groups are restraint to the value in qin
        ### self.noh_frozen     : Heavy atoms in the whole molecule are restraint to
        ###                       their value in qin
        ### self.h_groups_frozen: H atoms in groups are restraint to their value
        ###                       in the qin file.
        ### self.h_equiv        : Make degenerate hydrogen atoms equal
        ### self.all_equiv      : Make degenerate atoms (all) equal
        ### self.unfreeze_all   : Remove all restraints

        if self.stage == 'resp1_surr':

            _tmp_str.append(' nmol    = %d,\n' %self._mol_count)
            _tmp_str.append(' ihfree  = 1,\n')
            _tmp_str.append(' ioutopt = 1,\n')
            _tmp_str.append(' iqopt   = 1,\n')
            _tmp_str.append(' qwt     = 0.00500\n')

            self.groups_frozen   = False
            self.noh_frozen      = False
            self.h_groups_frozen = False
            self.h_equiv         = False
            self.all_equiv       = True
            self.unfreeze_all    = False

        elif self.stage == 'resp2_surr':

            _tmp_str.append(' nmol    = %d,\n' %self._mol_count)
            _tmp_str.append(' ihfree  = 1,\n')
            _tmp_str.append(' ioutopt = 1,\n')
            _tmp_str.append(' iqopt   = 2,\n')
            _tmp_str.append(' qwt     = 0.00100\n')

            self.groups_frozen   = True
            self.noh_frozen      = True
            self.h_groups_frozen = False
            self.h_equiv         = True
            self.all_equiv       = True
            self.unfreeze_all    = False

        elif self.stage == 'resp1':

            _tmp_str.append(' nmol    = %d,\n' %self._mol_count)
            _tmp_str.append(' ihfree  = 1,\n')
            _tmp_str.append(' ioutopt = 1,\n')
            _tmp_str.append(' iqopt   = 1,\n')
            _tmp_str.append(' qwt     = 0.00050\n')

            self.groups_frozen   = False
            self.noh_frozen      = False
            self.h_groups_frozen = False
            self.h_equiv         = False
            self.all_equiv       = True
            self.unfreeze_all    = False

        elif self.stage == 'resp2':

            _tmp_str.append(' nmol    = %d,\n' %self._mol_count)
            _tmp_str.append(' ihfree  = 1,\n')
            _tmp_str.append(' ioutopt = 1,\n')
            _tmp_str.append(' iqopt   = 2,\n')
            _tmp_str.append(' qwt     = 0.00100\n')

            self.groups_frozen   = False
            self.noh_frozen      = True
            self.h_groups_frozen = False
            self.h_equiv         = True
            self.all_equiv       = True
            self.unfreeze_all    = False

        _tmp_str.append('\n')
        _tmp_str.append(' &end')

        return ''.join(_tmp_str)

    def write(self, path):

        with open(path, 'w') as resp_file:

            resp_file.write(self.get_title())
            resp_file.write('\n')
            resp_file.write(self.get_cntrl())
            resp_file.write('\n')
            resp_file.write(self.get_mol())
            resp_file.write(self.get_group())
            resp_file.write(self.get_intermolecular())
            resp_file.write(self.get_atom_equiv())
            resp_file.write('\n')

        return True


def resp_surrogate(conn_surr,
    surr_conf_dict,
    surr_path_main,
    conn_i,
    stdout,
    stderr,
    cap_site=None,
    capsum=False):

    ### The input files for the resp program will be generated
    ### with the write_resp_in file generator.
    surr_resp_in = write_resp_in()

    ### Generate an esp file with all atomic and esp centers. Also we
    ### want to have the path to one representative mol2file from
    ### the qm calculation. This is mainly for getting the correct
    ### atom numbering.
    esp_list = list()

    mol_surr = conn_surr.get_surrogate_cap()

    surr_path_resp  = surr_path_main+"/resp"
    surr_path_esp   = surr_path_resp+"/esp"
    surr_path_resp1 = surr_path_resp+"/resp1"
    surr_path_resp2 = surr_path_resp+"/resp2"

    if not os.path.exists(surr_path_resp):
        os.mkdir(surr_path_resp)

    for conf_i in range(surr_conf_dict[conn_i]):
        conf_i_path = surr_path_main+"/conf%d/frag%d-conf%d_esp.log" %(conf_i, conn_i, conf_i)
        esp_path, mol2_path = preprocess_esp(conf_i_path)
        ### Match the mol2 structure obtained from the QM calculation
        ### with the rdkit mol object stored in the database. Perform
        ### some sanity checks on the match lists.
        mol_qm_surr  = Chem.MolFromMol2File(mol2_path, removeHs=False)
        mol_qm_surr  = fix_groups(mol_qm_surr)
        mol_surr     = fix_groups(mol_surr)
        matches      = mol_qm_surr.GetSubstructMatches(mol_surr)
        if len(matches) == 0:
            raise Exception("Did not any substructure matches.")
        if len(matches) > 1:
            raise Exception("Found %d substructure matches. Should be one only." %(len(matches)))
        if len(matches[0]) != mol_surr.GetNumAtoms():
            raise Exception("Found %d atoms in substructure matche. Expected %d." %(len(matches[0]),\
                        mol_surr.GetNumAtoms()))
        esp_list.append(esp_path)
        surr_resp_in.add_mol(mol_qm_surr, Chem.GetFormalCharge(mol_qm_surr), 1.0, "Conf %d" %conf_i)
        cat_esp(esp_list,surr_path_esp)
        ### Generate intermolecular restraints
        surr_resp_in.add_intermolecular(0,conf_i)

    if capsum:
        capsum_atom_i_list = list()
        capsum_mol_i_list  = list()

    matches = matches[0]
    ### Generate the group restraints
    atom_i_list = list()
    mol_i_list  = list()
    cap_neighbors = list()
    for cap_idx in conn_surr.get_rcap_map():
        cap_neighbors.append(matches[cap_idx])
    for cap_idx in conn_surr.get_rcap_map():
        atm_idx = matches[cap_idx]
        atm     = mol_qm_surr.GetAtomWithIdx(atm_idx)
        atom_i_list.append(atm_idx)
        mol_i_list.append(0)
        for neighbor in atm.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor.GetAtomicNum() == 1 \
            and neighbor_idx not in cap_neighbors:
                atom_i_list.append(neighbor_idx)
                mol_i_list.append(0)
    if not capsum:
        surr_resp_in.add_group(atom_i_list, mol_i_list, 0.0)
    else:
        capsum_atom_i_list += atom_i_list
        capsum_mol_i_list  += mol_i_list

    atom_i_list = list()
    mol_i_list  = list()
    cap_neighbors = list()
    for cap_idx in conn_surr.get_lcap_map():
        cap_neighbors.append(matches[cap_idx])
    for cap_idx in conn_surr.get_lcap_map():
        atm_idx = matches[cap_idx]
        atm     = mol_qm_surr.GetAtomWithIdx(atm_idx)
        atom_i_list.append(atm_idx)
        mol_i_list.append(0)
        for neighbor in atm.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor.GetAtomicNum() == 1 \
            and neighbor_idx not in cap_neighbors:
                atom_i_list.append(neighbor_idx)
                mol_i_list.append(0)
    if not capsum:
        surr_resp_in.add_group(atom_i_list, mol_i_list, 0.0)
    else:
        capsum_atom_i_list += atom_i_list
        capsum_mol_i_list  += mol_i_list

        surr_resp_in.add_group(capsum_atom_i_list, capsum_mol_i_list, 0.0)

    surr_resp_in.set_title('Stage 1 resp fit for connector id %d.' %conn_i)

    surr_resp_in.set_stage('resp1_surr')
    surr_resp_in.write(surr_path_resp1)

    surr_resp_in.set_title('Stage 2 resp fit. Surrogate cap connector id %d.' %conn_i)

    surr_resp_in.set_stage('resp2_surr')
    surr_resp_in.write(surr_path_resp2)

    resp_args = [resp_exe, "-O", "-i", surr_path_resp1, 
                           "-o", surr_path_resp+"/resp1.out",
                           "-e", surr_path_esp, 
                           "-t", surr_path_resp+"/resp1.crg", 
                           "-p", surr_path_resp+"/punch1",
                           "-s", surr_path_resp+"/esout1"]
    call(resp_args, stdout=stdout, stderr=stderr)

    resp_args = [resp_exe, "-O", "-i", surr_path_resp2, 
                           "-o", surr_path_resp+"/resp2.out",
                           "-q", surr_path_resp+"/resp1.crg",
                           "-e", surr_path_esp, 
                           "-t", surr_path_resp+"/resp2.crg", 
                           "-p", surr_path_resp+"/punch2",
                           "-s", surr_path_resp+"/esout2"]
    call(resp_args, stdout=stdout, stderr=stderr)

    ### Finally, write mol2 files for each conformer with the 
    ### calculated resp charges
    for conf_i in range(surr_conf_dict[conn_i]):
        mol2_path     = surr_path_main+"/conf%d/frag%d-conf%d_esp.mol2" \
                                       %(conf_i, conn_i, conf_i)
        mol2_out_path = surr_path_main+"/conf%d/frag%d-conf%d_resp.mol2"\
                                       %(conf_i, conn_i, conf_i)

        ante_args = [ante_exe, "-i", mol2_path, 
                               "-fi", "mol2", 
                               "-o", mol2_out_path, 
                               "-fo", "mol2",
                               "-cf", surr_path_resp+"/resp2.crg", 
                               "-c", "rc", 
                               "-at", "sybyl",
                               "-pf", "y",
                               "-dr", "no"]
        call(ante_args, stdout=stdout, stderr=stderr)

    return mol2_out_path, mol_qm_surr, matches