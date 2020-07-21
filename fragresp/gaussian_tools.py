import os

from subprocess import call
from rdkit import Chem

from fragresp.utils import logger
from fragresp.aux_progs import ante_exe

def get_freq(g09_log, freq_list):

    found_freq = False

    with open(g09_log, 'r') as f:
        for l in f:
            line = l.rstrip().split()
            if len(line) == 0:
                continue
            if line[0] == "Frequencies":
                freq_list.append(float(line[2]))
                freq_list.append(float(line[3]))
                freq_list.append(float(line[4]))

                found_freq = True

    return found_freq

def is_normal_terminate(g09_log):

    is_normal = False

    with open(g09_log, 'r') as f:
        for l in f:
            line = l.rstrip().split()
            if len(line) == 0:
                continue
            if len(line) > 4:
                if line[0] == "Normal" \
                and line[1] == "termination" \
                and line[2] == "of" \
                and line[3] == "Gaussian":
                    is_normal = True

    return is_normal

def check_opt(conf_dict,
              mol_list,
              qm_dir=".",
              logfile="check_qm.log",
              stdout=None,
              stderr=None):

    if stdout==None:
        stdout=open(os.devnull, 'w')
    if stderr==None:
        stderr=open(os.devnull, 'w')

    check_log = logger(logfile)
    check_log.log("# Logfile QM check")

    for frag_i, conf_N in conf_dict.items():

        check_log.log("Checking frag%s..." %frag_i)

        mainpath = qm_dir+"/"+"frag%d" %frag_i

        ante_args_general = ["-fi", "gout", "-pf", "y", "-at", "sybyl", "-dr", "no",]

        if not os.path.exists(mainpath):
            check_log.log(" ! Warning ! %s not found." %mainpath)
            continue

        for conf_i in range(conf_N):

            check_log.log("Checking frag%d-conf%d..." %(frag_i, conf_i))

            fragname = 'frag%d-conf%d' %(frag_i,conf_i)

            fragpath = mainpath+'/'+'conf%d' %conf_i
            logpath  = fragpath+"/"+fragname+"_opt.log"
            if not os.path.exists(fragpath):
                check_log.log(" ! Warning ! %s not found." %fragpath)
                continue

            if not os.path.exists(logpath):
                check_log.log(" ! Warning ! %s not found" %logpath)
                continue

            freq_list  = list()
            if get_freq(logpath, freq_list):
                freq_str   = 'Frequencies [cm**-1]: '
                freq_ok    = True
                for freq_i, freq in enumerate(freq_list):
                    if freq < 0.:
                        check_log.log(" ! Warning ! Imaginary frequency %f cm**-1 for mode %d" %(freq, freq_i))
                        freq_ok = False
                    freq_str += str(freq)
                    freq_str += ' '
                if freq_ok:
                    check_log.log(" Frequencies OK!")
                #check_log.log(freq_str)
            else:
                check_log.log(" ! Warning ! No frequencies found in %s" %logpath)
            del freq_list

            if is_normal_terminate(logpath):
                check_log.log(" Normal termination of Gaussian.")
            else:
                check_log.log(" ! Warning ! An error occured. Check Gaussian log.")

            mol2path  = fragpath+"/"+fragname+"_opt.mol2"
            ante_args = [ante_exe, "-i", logpath, "-o", mol2path, "-fo", "mol2"] + ante_args_general

            call(ante_args, stdout=stdout, stderr=stderr)

            mol_opt = Chem.MolFromMol2File(mol2path, removeHs=False)
            matches = mol_opt.GetSubstructMatches(mol_list[frag_i])
            if len(matches) > 1:
                check_log.log(" ! Warning ! Found more than one substructure match.")
            if len(matches) == 0:
                check_log.log(" ! Warning ! Did not find any substructure match.")

    check_log.close()

def get_esp(g09_log, esp_vals, esp_crds):

    found_esp_vals = False
    found_esp_crds = False
    first_esp_idx  = ""

    with open(g09_log, 'r') as f:
        for l in f:
            line = l.rstrip().split()
            if len(line) == 0:
                continue
            if len(line)>3:
                if line[0]  == "ESP" \
                and line[1] == "Fit" \
                and line[2] == "Center":
                    esp_crds.append(list())
                    esp_crds[-1].append(float(line[6]))
                    esp_crds[-1].append(float(line[7]))
                    esp_crds[-1].append(float(line[8]))
                    if not found_esp_crds:
                        first_esp_idx = line[3]
                    found_esp_crds = True

            if len(line)>2:
                if line[1] == "Fit" \
                and (found_esp_vals or line[0] == first_esp_idx):
                    esp_vals.append(float(line[2]))
                    found_esp_vals = True

    return found_esp_vals*found_esp_crds

def check_esp(conf_dict,
              qm_dir=".",
              logfile="check_esp.log"):

    check_log = logger(logfile)
    check_log.log("# Logfile QM ESP check")

    for frag_i, conf_N in conf_dict.items():

        check_log.log("Checking frag%s..." %frag_i)

        mainpath = qm_dir+"/"+"frag%d" %frag_i

        if not os.path.exists(mainpath):
            check_log.log(" ! Warning ! %s not found." %mainpath)
            continue

        for conf_i in range(conf_N):

            check_log.log("Checking frag%d-conf%d..." %(frag_i, conf_i))

            fragname = 'frag%d-conf%d' %(frag_i,conf_i)

            fragpath = mainpath+'/'+'conf%d' %conf_i
            logpath  = fragpath+"/"+fragname+"_esp.log"
            if not os.path.exists(fragpath):
                check_log.log(" ! Warning ! %s not found." %fragpath)
                continue

            if not os.path.exists(logpath):
                check_log.log(" ! Warning ! %s not found" %logpath)
                continue

            
            esp_vals_list = list()
            esp_crds_list = list()
            if get_esp(logpath, esp_vals_list, esp_crds_list):
                n_crds = len(esp_crds_list)
                n_vals = len(esp_vals_list)
                if n_crds == n_vals:
                    check_log.log(" Found %d esp fit centers." %n_vals)
                else:
                    check_log.log("! Warning ! Found %d fit centers but only %d fit values." %(n_crds, n_vals))

            if is_normal_terminate(logpath):
                check_log.log(" Normal termination of Gaussian.")
            else:
                check_log.log(" ! Warning ! An error occured. Check Gaussian log.")

    check_log.close()

def get_energy(g09_log, ene_list):

    found_ene = False
    with open(g09_log, 'r') as f:
        for l in f:
            line = l.rstrip().split()
            if len(line) < 2:
                continue
            if line[0] == "SCF" and line[1] == "Done:":
                ene_list.append(float(line[4]))
                found_ene = True

    return found_ene