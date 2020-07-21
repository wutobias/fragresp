import argparse
from fragresp.scripts.run_modes import run_modes

def parse_args():

    parser = argparse.ArgumentParser(
            description="Executable for FragResp Version 0.4. A processing pipeline for automated derivation of RESP charges \
within large sets of molecules. Written by Tobias Hüfner, tobias.wulsdorf@gmail.com")

    only_decompose ="Only used in --mode==decompose."
    only_conformers="Only used in --mode=conformers."
    only_check_opt ="Only used in --mode=check_opt."
    only_check_esp ="Only used in --mode=check_esp."
    only_resp      ="Only used in --mode=resp."
    only_datadump  ="Only used in --mode=datadump."

    no_default_arg ="There is not default value for this argument."
    no_default_opt ="There is not default value for this option."

    ### Required Arguments
    ### ~~~~~~~~~~~~~~~~~~

    required = parser.add_argument_group('required arguments')
    required.add_argument('-i', '--input', 
        required=True, 
        type=str, 
        default=None,
        help="This is the input file. It contains all the paramters that should remain unchanged during \
the process. %s" %no_default_opt)

    required.add_argument('-m', '--mode',
        required=True,
        type=str,
        default=None,
        help="With this option the current run mode is specified. Note, that each run mode has specific \
additional arguments. %s" %no_default_opt,
        choices=run_modes.keys())

    ### Optional Arguments ###
    ### ~~~~~~~~~~~~~~~~~~ ###

    parser._action_groups.append(parser._action_groups.pop(1))

    parser.add_argument('-l', '--linker', 
        required=False, 
        type=str, 
        default=None,
        help="This is the linker file. The content of this file instruct the decomposition routine how to \
break bonds and add capping groups to the created fragments. Also, in this file surrogate molecules are the mapping \
onto the capping groups are defined. %s %s" %(no_default_arg, only_decompose))

    parser.add_argument('-mol', '--mol2list', 
        required=False, 
        type=str, 
        default=None,
        help="This is the file that contains the (unique) names and paths to the molecules that are going to be \
decomposed. Note, that the molecules must be in mol2-format with sybyl atom types, other formats are currently not implemented but might \
come in future releases. It is possible to subsequently use different mol2lists together with the --add argument. %s %s" %(no_default_arg, only_decompose))

    parser.add_argument('-db', '--dumpdir',
        required=False,
        type=str,
        default="datadump",
        help='Path to datadump directory. %s' %only_datadump)

    parser.add_argument('-pl', '--pipeline',
        required=False,
        type=str,
        default="openeye",
        choices=["openeye", "openbabel"],
        help='Pipeline for conformer generation. Note, that the OpenEye pipeline requires a valid licence file, stored under $OE_LICENSE. Default is \
\'openeye\' (=use the OpenEye pipeline). %s' %only_conformers)

    parser.add_argument('-a', '--add',
        action='store_true',
        help='Add molecules to existing database. Default is off (=do not add molecules to the database). %s' %only_decompose)

    parser.add_argument('-lim', '--limit',
        required=False,
        type=int,
        default=3,
        help='Sets the number of diverse conformers to select from the least interacting conformers set in the second step of the ELF algorithm. \
%s' %only_conformers)

    parser.add_argument('-per', '--percentage',
        required=False,
        type=float,
        default=50.,
        help='Sets the percentage of input conformers to select as the electrostatically least-interacting set in the first step of the ELF algorithm. \
%s' %only_conformers)

    parser.add_argument('-fl', '--fraglog',
        required=False,
        type=str,
        default="frag_prep.log",
        help='Log file for conformer generation and QM preparation for fragments. %s' %only_conformers)

    parser.add_argument('-sl', '--surrlog',
        required=False,
        type=str,
        default="surr_prep.log",
        help='Log file for conformer generation and QM preparation for surrogates. %s' %only_conformers)

    parser.add_argument('-ml', '--mollog',
        required=False,
        type=str,
        default="mol_prep.log",
        help='Log file for conformer generation and QM preparation for molecules without \
decomposition. %s' %only_conformers)

    parser.add_argument('-q', '--queue',
        required=False,
        type=str,
        default="marc2",
        choices=['marc2', 'slurm', 'condor', 'none'],
        help='Queue used for gaussian09 jobs. Default is marc2. %s' %only_conformers)

    parser.add_argument('-fs', '--fittosurr',
        action='store_true',
        help='Derives the capping charges from the surrogates. Note, that the capping charges then are not allowed to vary freely. \
Default is off (=do not fit to surrogate). %s' %only_resp)

    parser.add_argument('-by', '--bycharge',
        action='store_true',
        help='Perform seperate resp charge calculation for fragments with different total charge. The default for this \
option is off (=do not perform seperate resp calculations for fragments with different total charge). %s' %only_resp)

    parser.add_argument('-cs', '--capsum',
        action='store_true',
        help='If this is activated, the charges of terminal linker groups are constraint, such that the charges of *both* capping groups sum up to zero. \
If this option is deactivated, then the charges of the rcap group are used directly for the linker. The default for this \
option is off (=do not constraint the charges of both capping groups to sum up to zero. Use the charges of the rcap group directly.)%s' %only_resp)

    parser.add_argument('-re', '--remap',
        action='store_true',
        help='Maps the derived charges back to the original set of molecules. Note, that if this option is activated, then the molecules \
will be stored as mol2 files in --remap_dir. The default for this option is off (=do not remap to original molecules). %s' %only_resp)

    parser.add_argument('-rd', '--remapdir',
        required=False,
        type=str,
        default="remap_dir",
        help='Path for storing mol2 files after remapping fragment charges. %s' %only_resp)

    parser.add_argument('-v', '--verbose',
        action='store_true',
        help='Verbosity output. Note, that this generates *a lot* of output in --mode=decompose. Default is off.')

    parser.add_argument('-d', '--debug',
        action='store_true',
        help='Run in debug mode. Default is off.')

    args = parser.parse_args()

    if args.mode == 'decompose':
        if args.mol2list == None:
            parser.error("In --mode=decompose, '--mol2list' cannot be empty.")
        if args.linker == None and not args.add:
            parser.error("In --mode=decompose without '--add' argument, '--linker' cannot be empty.")

    return args