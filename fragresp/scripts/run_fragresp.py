import os
from fragresp.scripts.parse_args import parse_args
from fragresp.scripts.parse_args import run_modes
from fragresp.scripts.utils import read_infile
from fragresp.scripts.utils import read_initfile
from fragresp.scripts.utils import parse_linker

def main():

    options_infile = list()
    options_infile.append('frag_dir')
    options_infile.append('mol_dir')
    options_infile.append('surr_cap_dir')
    options_infile.append('database')

    options_linker = list()
    options_linker.append('name')
    options_linker.append('linker')
    options_linker.append('lcap')
    options_linker.append('lanc')
    options_linker.append('lancmap')
    options_linker.append('rcap')
    options_linker.append('ranc')
    options_linker.append('rancmap')
    options_linker.append('breakring')
    options_linker.append('surrogate')
    options_linker.append('s_lanc')
    options_linker.append('s_ranc')
    options_linker.append('terminal')
    options_linker.append('terminal_anc')

    #Must be empty!
    options_mol2 = list()

    args           = parse_args()
    mol2_list      = None
    linker_objects = None

    args_infile  = read_infile(args.input, options_infile)
    args_infile['database']     = args_infile['database'][0]
    args_infile['frag_dir']     = args_infile['frag_dir'][0]
    args_infile['mol_dir']      = args_infile['mol_dir'][0]
    args_infile['surr_cap_dir'] = args_infile['surr_cap_dir'][0]

    if args.debug:
        print ("--args")
        for key,item in args_infile.items():
            print ("  %s: %s" %(key,item))

    if not args.add and args.mode=='decompose':
        linker_list = read_initfile(args.linker, options_linker)
        if args.debug:
            print ("--linkerlist")
            for i, linker in enumerate(linker_list):
                print ("  linker %d" %i)
                for key,item in linker.items():
                    print ("  %s: %s" %(key,item))
        linker_objects = parse_linker(linker_list, verbose=args.verbose)
    

    if args.mode=='decompose':
        mol2_list   = read_initfile(args.mol2list, options_mol2)
        if args.debug:
            print ("--mol2list")
            for i, mol2 in enumerate(mol2_list):
                print ("  mol2list %d" %i)
                for key,item in mol2.items():
                    print ("  %s: %s" %(key,item))

    if args.mode == 'decompose':
        run_modes['decompose'](mol2_list, 
            args_infile['database'], 
            linker_objects, 
            args.add, 
            args.verbose)

    if args.mode == 'datadump':
        run_modes['datadump'](args_infile['database'], 
            args.dumpdir)

    if args.mode == 'conformers':
        run_modes['conformers'](args_infile['database'],
            args.add,
            args.limit,
            args.percentage,
            args_infile['frag_dir'],
            args_infile['surr_cap_dir'],
            args_infile['mol_dir'],
            args.fraglog,
            args.surrlog,
            args.mollog,
            args.pipeline,
            args.queue)

    if args.mode == 'check_opt':
        stdout=open(os.devnull, 'w')
        stderr=open(os.devnull, 'w')
        if args.verbose or args.debug:
            stdout=open("stdout", "w")
            stderr=open("stderr", "w")
        run_modes['check_opt'](args_infile['database'],
            args_infile['frag_dir'],
            args_infile['surr_cap_dir'],
            args_infile['mol_dir'],
            stdout,
            stderr)

    if args.mode == 'check_esp':
        run_modes['check_esp'](args_infile['database'],
            args_infile['frag_dir'],
            args_infile['surr_cap_dir'],
            args_infile['mol_dir'])

    if args.mode == 'resp':
        stdout=open(os.devnull, 'w')
        stderr=open(os.devnull, 'w')
        if args.verbose or args.debug:
            stdout=open("stdout", "w")
            stderr=open("stderr", "w")
        run_modes['resp'](args_infile['database'],
            args_infile['frag_dir'],
            args_infile['surr_cap_dir'],
            args_infile['mol_dir'],
            args.fittosurr,
            args.remap,
            args.bycharge,
            args.capsum,
            args.remapdir,
            args.energy_weighting,  
            stdout,
            stderr)

def entry_point():

    main()

if __name__ == '__main__':

    entry_point()
