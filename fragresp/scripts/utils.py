import copy
from collections import OrderedDict
from fragresp.data_strc import connector

def read_infile(path, option_choices):

    option_dict    = OrderedDict()

    with open(path, 'r') as file:

        for line in file:
            l = line.rstrip().lstrip().split()
            
            if len(l) == 0:
                continue
            if l[0].startswith('#'):
                continue

            if len(l) == 1:
                raise IOError("Option %s not understood." %l[0])
            elif l[0] not in option_choices:
                raise IOError("Option %s not known." %l[0])
            option_dict[l[0]] = list()
            for val in l[1:]:
                option_dict[l[0]].append(val)

    return option_dict


def read_initfile(path, option_choices):

    init_list    = list()

    init=False

    with open(path, 'r') as file:

        for line in file:
            l = line.rstrip().lstrip().split()
            
            if len(l) == 0:
                continue
            if l[0].startswith('#'):
                continue

            if l[0]=='init' and init:
                raise IOError("Linker section must be closed before new section can be started.")

            if l[0]=='end' and not init:
                raise IOError("Must open linker section before start a new one.")

            if l[0]=='init':
                init=True
                init_list.append(OrderedDict())
                continue

            if l[0]=='end':
                init=False
                continue

            if init:
                if len(l) == 1:
                    raise IOError("Option %s not understood." %l[0])
                elif len(option_choices)>0 and l[0] not in option_choices:
                    raise IOError("Option %s not known." %l[0])
                init_list[-1][l[0]] = list()
                for val in l[1:]:
                    init_list[-1][l[0]].append(val)

    return init_list

def parse_linker(linker_dict_list, verbose=False):

    linker_object_list = list()

    for linker_dict in linker_dict_list:

        if verbose:
            print ("Reading linker %s" %linker_dict['name'])

        tointlist = ['lanc', 'ranc', 'lancmap', 'rancmap', 's_lanc', 's_ranc']

        ring_save = True
        if linker_dict['breakring'][0] == 'on':
            ring_save = False

        terminal_linker = False
        if linker_dict['terminal'][0] == 'on':
            terminal_linker = True

        for anc in tointlist:
            if not (anc.startswith("s_") and terminal_linker):
                for i,n in enumerate(linker_dict[anc]):
                    linker_dict[anc][i] = int(n)
            elif terminal_linker:
                for i,n in enumerate(linker_dict["terminal_anc"]):
                    linker_dict["terminal_anc"][i] = int(n)

        conn = connector("".join(linker_dict['name']))
        conn.init_connector("".join(linker_dict['linker']))
        if not terminal_linker:
            conn.init_lcap("".join(linker_dict['lcap']))
            conn.init_rcap("".join(linker_dict['rcap']))
        conn.init_lanc(linker_dict['lanc'])
        conn.init_ranc(linker_dict['ranc'])
        conn.init_lancmap(linker_dict['lancmap'])
        conn.init_rancmap(linker_dict['rancmap'])
        conn.init_ring(ring_save)
        conn.init_terminal(terminal_linker)
        if terminal_linker:
            if 'terminal_anc' not in linker_dict:
                raise IOError("If terminal=on, the terminal_anc must be given.")
            conn.init_terminal_anc(linker_dict['terminal_anc'])
        else:
            conn.init_surrogate("".join(linker_dict['surrogate']))
            conn.init_s_lanc(linker_dict['s_lanc'])
            conn.init_s_ranc(linker_dict['s_ranc'])
            conn.init_surrogate_cap()
        if not conn.check(verbose):
            raise Warning("Linker %s not fully defined." %linker_dict['name'])

        linker_object_list.append(copy.copy(conn))

    return linker_object_list

example_input="""
database frag_db.pickle
frag_dir frag_dir
surr_cap_dir surr_cap_dir
"""

example_mol2list="""
### mol2 file list ###

init mol2
1K21 mol2/1K21.mol2
2ZC9 mol2/2ZC9.mol2
2ZDA mol2/2ZDA.mol2
2ZDV mol2/2ZDV.mol2
2ZF0 mol2/2ZF0.mol2
2ZFF mol2/2ZFF.mol2
2ZFP mol2/2ZFP.mol2
2ZGX mol2/2ZGX.mol2
2ZO3 mol2/2ZO3.mol2
3BIU mol2/3BIU.mol2
3BIV mol2/3BIV.mol2
3DHK mol2/3DHK.mol2
3DUX mol2/3DUX.mol2
end mol2

"""

example_linker="""
### Linker file for use with Thrombin inhibitors ###

init
name Primary Peptide
linker CNC(=O)C
# NME cap
lcap NC
lanc 2 3 4
lancmap 1 2 -1
# ACE cap
rcap C(=O)C
ranc 0 1
rancmap -1 5
breakring off
# Alanine surrogate
surrogate NC(=O)[C@@H](C)NC(=O)
s_lanc 0
s_ranc 6 7
end

init
name  Secondary Peptide
linker CN(C)C(=O)C
# NME2 cap
lcap N(C)C
lanc 3 4 5
lancmap 1 2 -1
# ACE
rcap C(=O)C
ranc 0 1 2
rancmap -1 5 -1
breakring off
# Me-Alanine surrogate
surrogate NC(=O)[C@@H](C)N(C)C(=O)
s_lanc 0
s_ranc 7 8
end

init
name Primary Sulfonamide
linker CNS(=O)(=O)C
# Sulfo-NME cap
lcap NC
lanc 2 3 4 5
lancmap 1 2 -1 -1
# Sulfo-ACE cap
rcap S(=O)(=O)C
ranc 0 1
rancmap -1 6
breakring off
# Sulfon amide Alanine
surrogate NS(=O)(=O)[C@@H](C)NS(=O)(=O)
s_lanc 0
s_ranc 7 8 9
end

"""