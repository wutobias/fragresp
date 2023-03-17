import os
from openeye.oechem import *
from openeye.oeomega import *
from openeye.oeff import *
from openeye.oequacpac import *
from openeye.oeszybki import *

from rdkit import Chem
from fragresp.utils import logger

def make_g09(name,
             conf,
             nproc=4,
             mem=1900,
             queue='marc2'):

    queue_choices = ['marc2', 'slurm', 'condor', 'none']
    if queue not in queue_choices:
        print ("queue %s not known. Using queue 'none'.")
        queue = 'none'

    radius_dict = {'I'  : 1.98}

    sumFCharge    = 0
    com_crds      = ''
    com_opt       = ''
    com_esp       = ''
    com_basis     = ''
    crds          = conf.GetCoords()

    lelement_list     = list()
    helement_list     = list()
    has_user_basisset = False
    for atm in conf.GetAtoms():
        atm_number  = atm.GetAtomicNum()
        element     = OEGetAtomicSymbol(atm_number)
        if atm_number > 35:
            if not element in helement_list:
                helement_list.append(element)
            has_user_basisset = True
        else:
            if not element in lelement_list:
                lelement_list.append(element)
        sumFCharge += atm.GetFormalCharge()
        com_crds   += '   '
        com_crds   += element
        com_crds   += '   '
        for crd in crds[atm.GetIdx()]:
            if crd > 0.:
                com_crds += " "
            com_crds += '%5.10f ' %crd
        com_crds += '\n'
    
    if queue in ['none', 'slurm']:
        com_opt += '%%nproc=%d\n'   %nproc
        com_opt += '%%mem=%dMB\n'   %mem
    com_opt += '%%chk=%s.chk\n' %name
    if has_user_basisset:
        com_opt += '#B3LYP/GenECP 5d 7f pseudo=read '
    else:
        com_opt += '#B3LYP/6-31G* '
    com_opt += 'Integral=(Grid=UltraFine) freq Opt pop=none\n'
    com_opt += '\n'
    com_opt += 'Optimization for %s\n' %name
    com_opt += '\n'
    com_opt += '%d  1\n' %int(sumFCharge)
    com_opt += com_crds
    com_opt += '\n'
    if has_user_basisset:
        for element in helement_list:
            com_basis += '%s ' %element
        com_basis += '0\n'
        com_basis += 'sdd\n'
        com_basis += '****\n'
        for element in lelement_list:
            com_basis += '%s ' %element
        com_basis += '0\n'
        com_basis += '6-31G*\n'
        com_basis += '****\n'
        com_basis += '\n'
        for element in helement_list:
            com_basis += '%s ' %element
        com_basis += '0\n'
        com_basis += 'sdd\n'
        com_basis += '\n'
    com_opt += com_basis
    com_opt += '\n'

    if queue in ['none', 'slurm']:
        com_esp += '%%nproc=%d\n'  %nproc
        com_esp += '%%mem=%dMB\n'  %mem
    com_esp += '%%chk=%s.chk\n'    %name
    if has_user_basisset:
        com_esp += '#HF/GenECP 5d 7f pseudo=read '
    else:
        com_esp += '#HF/6-31G* '
    com_esp += 'Geom=AllCheck Guess=Read Pop=(MK,ReadRadii) '
    com_esp += 'Integral=(Grid=UltraFine) iop(6/33=2)\n'
    com_esp += '\n'
    com_esp += com_basis
    if has_user_basisset:
        for element in helement_list:
            com_esp += '%s %f\n' %(element, radius_dict[element])
    com_esp += '\n'

    return com_opt, com_esp

def filter_conformers(mol, limit=3, percentage=50.):

    charges = OEAM1BCCELF10Charges()
    charges.SetLimit(limit)
    charges.SetPercentage(percentage)
    charges.SetReturnSelectedConfs(True)
    OEAssignCharges(mol, charges)


def generate_conformers(mol, optimize=True):

    omegaOpts = OEOmegaOptions()
    omega = OEOmega(omegaOpts)
    omega(mol)

    if optimize:
        optimize_conformers(mol)
    

def optimize_conformers(mol):

    opts = OESzybkiOptions()
    opts.GetOptOptions().SetOptimizerType(OEOptType_NEWTON)
    opts.GetGeneralOptions().SetForceFieldType(OEForceFieldType_MMFF94S)
    opts.GetSolventOptions().SetSolventModel(OESolventModel_Sheffield)
    opts.GetSolventOptions().SetChargeEngine(OEChargeEngineNoOp())

    sz  = OESzybki(opts)
    res = OESzybkiResults()

    for conf in mol.GetConfs():
        sz(conf, res)

def prep_qm(frag_list,
            overwrite,
            limit,
            percentage,
            write_mol2=True,
            nproc=2,
            mem=1900,
            queue='none',
            qm_dir=".",
            opt_batch="submit_opt_g09.sh",
            esp_batch="submit_esp_g09.sh",
            logfile="prep_qm.log",
            include_list=None):

    conf_list = list()

    if include_list==None:
        include_list=range(len(frag_list))

    prep_log = logger(logfile)
    prep_log.log("### OpenEye conformer pipeline")

    if not os.path.exists(qm_dir):
        os.mkdir(qm_dir)

    if queue:
        g09_cmd = 'g09'
    else:
        from pkg_resources import resource_filename
        g09_cmd = resource_filename("fragresp.data", f"submit-scripts/{queue}.py")

    opt_batch_file = logger(qm_dir+"/"+opt_batch)
    esp_batch_file = logger(qm_dir+"/"+esp_batch)

    opt_batch_file.log("#!/bin/bash")
    opt_batch_file.log()
    esp_batch_file.log("#!/bin/bash")
    esp_batch_file.log()
    for frag_i in include_list:

        frag = frag_list[frag_i]

        mainpath  = qm_dir+"/"+"frag%d" %frag_i

        if not os.path.exists(mainpath):
            os.mkdir(mainpath)
        elif not overwrite:
            conf_list.append(-1)
            continue

        ### Always isomericSmiles=True, otherwise omega does not
        ### know what to do.
        smi = Chem.MolToSmiles(frag, isomericSmiles=True)

        prep_log.log("Fragment %d" %frag_i)
        prep_log.log("Smi      %s" %smi)

        mol = OEGraphMol()
        OESmilesToMol(mol, smi)
        mol = OEMol(mol)

        generate_conformers(mol)
        prep_log.log("initial conformers: %d" %mol.NumConfs())
        filter_conformers(mol, limit, percentage)
        prep_log.log("filtered conformers:%d" %mol.NumConfs())
        conf_list.append(mol.NumConfs())

        for conf_i, conf in enumerate(mol.GetConfs()):

            fragname = 'frag%d-conf%d' %(frag_i,conf_i)

            fragpath = mainpath+'/'+'conf%d' %conf_i
            subpath  = "frag%d" %frag_i + '/'+'conf%d' %conf_i
            if not os.path.exists(fragpath):
                os.mkdir(fragpath)        

            com_opt, com_esp = make_g09(fragname, conf, nproc, mem, queue)

            f = logger(fragpath+"/"+fragname+"_opt.com")
            f.log(com_opt)
            f.close()
            opt_batch_file.log("%s %s" %(g09_cmd, subpath+"/"+fragname+"_opt.com"))
            del f

            f = logger(fragpath+"/"+fragname+"_esp.com")
            f.log(com_esp)
            f.close()
            esp_batch_file.log("%s %s" %(g09_cmd, subpath+"/"+fragname+"_esp.com"))
            del f

        if write_mol2:
            ofs = oemolostream()
            ofs.open(mainpath+"/"+"f%d.mol2" %frag_i)
            OEWriteMolecule(ofs, mol)

            del ofs

        prep_log.log()

    prep_log.close()

    opt_batch_file.close()
    esp_batch_file.close()

    del opt_batch_file
    del esp_batch_file

    return conf_list