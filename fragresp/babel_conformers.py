import os
from openbabel import openbabel as ob
from openbabel import pybel as pb

from rdkit import Chem
from fragresp.utils import logger

def make_g09(name,
             conf,
             nproc=4,
             mem=1900,
             queue='marc2'):

    queue_choices = ['marc2', 'slurm', 'condor', 'none']
    if queue not in queue_choices:
        print( "queue %s not known. Using queue 'none'.")
        queue = 'none'

    radius_dict = {'I'  : 1.98}

    sumFCharge    = 0
    com_crds      = ''
    com_opt       = ''
    com_esp       = ''
    com_basis     = ''

    lelement_list     = list()
    helement_list     = list()
    has_user_basisset = False

    for atm in ob.OBMolAtomIter(conf):
        atm_number  = atm.GetAtomicNum()
        element     = ob.GetSymbol(atm_number)
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
        crds        = [atm.GetX(), atm.GetY(), atm.GetZ()]
        for crd in crds:
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
        com_opt += '%%nproc=%d\n'   %nproc
        com_opt += '%%mem=%dMB\n'   %mem
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
    
    ### This is an (approximate) reimplementation of the ELF
    ### algorithm as implemented in the OpenEye Toolkits.
    
    ene_init = list()
    ene_elf  = list()
    
    percentage /= 100.
    if percentage>1.:
        percentage=1.
    elif percentage<0.:
        percentage=0.1
    
    n_confs = mol.NumConformers()
    mff     = ob.OBForceField.FindForceField("mmff94")
    for conf_i in range(n_confs):
        mol.SetConformer(conf_i)
        mff.Setup(mol)
        mff.GetCoordinates(mol)
        ene_init.append(mff.Energy())
    confs_init = sorted(range(n_confs), key=ene_init.__getitem__)
    cut_1st    = int(n_confs*percentage)
    if cut_1st==0:
        cut_1st=1
    if cut_1st<limit:
        limit=cut_1st
        
    for a in ob.OBMolAtomIter(mol):
        chg = a.GetPartialCharge()
        if chg < 0.:
            a.SetPartialCharge(chg*-1.)
    
    sff = ob.OBForceField.FindForceField("mmff94")
    for conf_i in range(cut_1st):
        mol.SetConformer(confs_init[conf_i])
        sff.Setup(mol)
        sff.GetCoordinates(mol)
        ene_elf.append(sff.Energy())
        
    confs_elf = sorted(range(cut_1st), key=ene_elf.__getitem__)
    confs_final = list()
    for conf_i in range(limit):
        confs_final.append(confs_init[confs_elf[conf_i]])
    del_i = 0
    for conf_i in range(n_confs):
        if conf_i not in confs_final:
            mol.DeleteConformer(del_i)
        else:
            del_i += 1


def generate_conformers(mol, optimize=True, conf_search='diverse'):

    if conf_search == 'diverse':
        rmsd_cutoff    = 1.0
        conf_cutoff    = 40000000
        energy_cutoff  = 10.0
        confab_verbose = False

        # Run Confab conformer generation
        ### see also:
        ### https://gist.github.com/kylebarlow/1756ea399ba6bfee3c2d3d054c17c3a3
        cff = ob.OBForceField.FindForceField("mmff94")
        cff.Setup(mol)
        cff.DiverseConfGen(rmsd_cutoff, conf_cutoff, energy_cutoff, confab_verbose)
        cff.GetConformers(mol)

    else:
        cs = ob.OBConformerSearch()
        numConformers = 30
        numChildren = 5
        mutability = 5
        convergence = 25
        cs.Setup(mol,numConformers,numChildren,mutability,convergence)
        cs.GetConformers(mol)
        
    if optimize:
        optimize_conformers(mol)


def optimize_conformers(mol):
        
    n_confs  = mol.NumConformers()
    mff      = ob.OBForceField.FindForceField("mmff94")
    for conf_i in range(n_confs):
        mol.SetConformer(conf_i)
        mff.Setup(mol)
        mff.SteepestDescent(100)
        mff.GetCoordinates(mol)


def prep_qm(frag_list,
            overwrite,
            limit,
            percentage,
            write_mol2=True,
            nproc=2,
            mem=1900,
            queue='marc2',
            qm_dir=".",
            opt_batch="submit_opt_g09.sh",
            esp_batch="submit_esp_g09.sh",
            logfile="prep_qm.log",
            include_list=None):

    conf_list = list()

    if include_list==None:
        include_list=range(len(frag_list))

    prep_log = logger(logfile)
    prep_log.log("### Openbabel conformer pipeline")

    if not os.path.exists(qm_dir):
        os.mkdir(qm_dir)

    if queue == 'marc2':
        g09_cmd = "subg09 -p %d -m %d -cc" %(nproc, mem)
    elif queue == 'slurm':
        g09_cmd = 'subg09_slurm.sh %d' %nproc
    elif queue == 'condor':
        g09_cmd = 'subg09_condor.sh 1'
    else:
        g09_cmd = 'g09'

    opt_batch_file = logger(qm_dir+"/"+opt_batch)
    esp_batch_file = logger(qm_dir+"/"+esp_batch)

    opt_batch_file.log("#!/bin/bash")
    opt_batch_file.log()
    opt_batch_file.log("_pwd=$PWD")
    opt_batch_file.log()
    esp_batch_file.log("#!/bin/bash")
    esp_batch_file.log()
    esp_batch_file.log("_pwd=$PWD")
    esp_batch_file.log()

    obConversion = ob.OBConversion()
    obConversion.SetInFormat("smi")
    obConversion.SetOutFormat("mol2")

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
        
        pbmol = pb.readstring("smi", smi)
        pbmol.make3D()
        mol   = pbmol.OBMol

        generate_conformers(mol)
        prep_log.log("initial conformers: %d" %mol.NumConformers())
        filter_conformers(mol, limit, percentage)
        prep_log.log("filtered conformers:%d" %mol.NumConformers())
        conf_list.append(mol.NumConformers())

        for conf_i in range(mol.NumConformers()):

            mol.SetConformer(conf_i)

            fragname = 'frag%d-conf%d' %(frag_i,conf_i)

            fragpath = mainpath+'/'+'conf%d' %conf_i
            subpath  = "frag%d" %frag_i + '/'+'conf%d' %conf_i
            if not os.path.exists(fragpath):
                os.mkdir(fragpath)        

            com_opt, com_esp = make_g09(fragname, mol, nproc, mem, queue)

            f = logger(fragpath+"/"+fragname+"_opt.com")
            f.log(com_opt)
            f.close()
            opt_batch_file.log("cd %s" %subpath)
            opt_batch_file.log("%s %s" %(g09_cmd, fragname+"_opt.com"))
            opt_batch_file.log("cd $_pwd")
            del f

            f = logger(fragpath+"/"+fragname+"_esp.com")
            f.log(com_esp)
            f.close()
            esp_batch_file.log("cd %s" %subpath)
            esp_batch_file.log("%s %s" %(g09_cmd, fragname+"_esp.com"))
            esp_batch_file.log("cd $_pwd")
            del f

        if write_mol2:
            obConversion.WriteFile(mol, mainpath+"/"+"f%d.mol2" %frag_i)

        prep_log.log()

    prep_log.close()

    opt_batch_file.close()
    esp_batch_file.close()

    del opt_batch_file
    del esp_batch_file

    return conf_list