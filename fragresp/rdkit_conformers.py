import os

import numpy as np

from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign
from fragresp.utils import logger
from fragresp.constants import max_confs

from pkg_resources import resource_filename

def make_g09(name,
             rdmol,
             conf_i,
             nproc=8,
             mem=1900):

    radius_dict = {'I'  : 1.98}

    sumFCharge    = 0
    com_crds      = ''
    com_opt       = ''
    com_esp       = ''
    com_basis     = ''

    lelement_list     = list()
    helement_list     = list()
    has_user_basisset = False

    conformer = rdmol.GetConformer(conf_i)
    positions = conformer.GetPositions()
    for atm in rdmol.GetAtoms():
        atm_number  = atm.GetAtomicNum()
        element     = atm.GetSymbol()
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
        crds        = positions[atm.GetIdx()].tolist()
        for crd in crds:
            if crd > 0.:
                com_crds += " "
            com_crds += '%5.10f ' %crd
        com_crds += '\n'

    com_opt += '%%nproc=%d\n'   %nproc
    com_opt += '%%mem=%dMB\n'   %mem
    com_opt += '%%chk=%s.chk\n' %name
    if has_user_basisset:
        com_opt += '#B3LYP/GenECP 5d 7f pseudo=read '
    else:
        com_opt += '#B3LYP/6-31G* '
    com_opt += 'Integral=(Grid=UltraFine) freq Opt=verytight pop=none\n'
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

    com_esp += '%%nproc=%d\n'   %nproc
    com_esp += '%%mem=%dMB\n'   %mem
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


def filter_conformers(rdmol, limit=3, percentage=50.):
    
    ### This is an (approximate) reimplementation of the ELF
    ### algorithm as implemented in the OpenEye Toolkits.

    from scipy.spatial import distance
    
    ene_init = list()
    ene_elf  = list()
    
    percentage /= 100.
    if percentage>1.:
        percentage=1.
    elif percentage<0.:
        percentage=0.1
    
    n_confs = rdmol.GetNumConformers()
    mp = AllChem.MMFFGetMoleculeProperties(rdmol, mmffVariant='MMFF94s')
    for conf_i in range(n_confs):
        mff = AllChem.MMFFGetMoleculeForceField(rdmol, mp, confId=conf_i)
        ene_init.append(mff.CalcEnergy())

    confs_init = sorted(range(n_confs), key=ene_init.__getitem__)
    cut_1st    = int(n_confs*percentage)
    if cut_1st==0:
        cut_1st=1
    if cut_1st<limit:
        limit=cut_1st

    charge_product = np.ones(
        (
            rdmol.GetNumAtoms(),
            rdmol.GetNumAtoms()
            ),
        dtype=float,
        )
    for atom_i in range(rdmol.GetNumAtoms()):
        q = mp.GetMMFFPartialCharge(atom_i)
        charge_product[atom_i:] *= q
        charge_product[:atom_i] *= q

    
    for conf_i in range(n_confs):
        conformer = rdmol.GetConformer(conf_i)
        positions = conformer.GetPositions()
        distance_matrix = distance.cdist(positions, positions)
        idxs = np.tril_indices(rdmol.GetNumAtoms(), -1)
        ene  = np.sum(
            1./distance_matrix[idxs] * charge_product[idxs]
            )
        ene_elf.append(ene)

    confs_elf = sorted(range(n_confs), key=ene_elf.__getitem__)
    to_delete = sorted(confs_elf[limit:], reverse=True)
    for conf_i in to_delete:
        rdmol.RemoveConformer(conf_i)


def generate_conformers(rdmol, optimize=True):

    param = rdDistGeom.ETKDGv2()
    param.pruneRmsThresh = 0.5
    cids = rdDistGeom.EmbedMultipleConfs(
        rdmol, 
        max_confs, 
        param)

    if optimize:
        AllChem.MMFFOptimizeMoleculeConfs(
            rdmol, 
            numThreads=0, 
            mmffVariant='MMFF94s'
            )


def optimize_conformers(rdmol):
        
    AllChem.MMFFOptimizeMoleculeConfs(
        rdmol, 
        numThreads=0, 
        mmffVariant='MMFF94s'
        )


def prep_qm(frag_list,
            overwrite,
            limit,
            percentage,
            write_sdf=True,
            nproc=8,
            mem=1900,
            queue='none',
            qm_dir=".",
            opt_batch="submit_opt_g09.sh",
            esp_batch="submit_esp_g09.sh",
            psi4_batch="submit_psi4.sh",
            logfile="prep_qm.log",
            include_list=None):

    conf_list = list()

    if include_list==None:
        include_list=range(len(frag_list))

    prep_log = logger(logfile)
    prep_log.log("### RDKit conformer pipeline")

    if not os.path.exists(qm_dir):
        os.mkdir(qm_dir)

    if queue == "none":
        g09_cmd = 'g09'
        psi4_cmd = 'run_psi4'
    else:
        from pkg_resources import resource_filename
        g09_cmd = resource_filename("fragresp.data", f"submit-gaussian/{queue}.py")
        psi4_cmd = resource_filename("fragresp.data", f"submit-psi4/{queue}.py")

    psi4_batch_file = logger(qm_dir+"/"+psi4_batch)
    opt_batch_file  = logger(qm_dir+"/"+opt_batch)
    esp_batch_file  = logger(qm_dir+"/"+esp_batch)

    psi4_batch_file.log("#!/bin/bash")
    psi4_batch_file.log()
    opt_batch_file.log("#!/bin/bash")
    opt_batch_file.log()
    esp_batch_file.log("#!/bin/bash")
    esp_batch_file.log()

    for frag_i in include_list:

        frag = frag_list[frag_i]
        Chem.SanitizeMol(frag)
        frag = AllChem.AddHs(frag, addCoords=True)

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
        
        generate_conformers(frag)
        prep_log.log("initial conformers: %d" %frag.GetNumConformers())
        filter_conformers(frag, limit, percentage)
        prep_log.log("filtered conformers:%d" %frag.GetNumConformers())
        conf_list.append(frag.GetNumConformers())

        for conf_i in range(frag.GetNumConformers()):

            fragname = 'frag%d-conf%d' %(frag_i,conf_i)

            fragpath = mainpath+'/'+'conf%d' %conf_i
            subpath  = "frag%d" %frag_i + '/'+'conf%d' %conf_i
            if not os.path.exists(fragpath):
                os.mkdir(fragpath)        

            com_opt, com_esp = make_g09(fragname, frag, conf_i, nproc, mem)

            f = logger(fragpath+"/"+fragname+"_opt.com")
            f.log(com_opt)
            f.close()
            opt_batch_file.log("%s %s %d" %(g09_cmd, subpath+"/"+fragname+"_opt.com", nproc))
            del f

            f = logger(fragpath+"/"+fragname+"_esp.com")
            f.log(com_esp)
            f.close()
            esp_batch_file.log("%s %s %d" %(g09_cmd, subpath+"/"+fragname+"_esp.com", nproc))
            del f

            psi4_batch_file.log("%s %s %d" %(psi4_cmd, subpath+"/"+fragname+"_opt.com", nproc))

        if write_sdf:
            w = Chem.SDWriter(mainpath+"/"+"f%d.sdf" %frag_i)
            res = list()
            mp = AllChem.MMFFGetMoleculeProperties(
                frag, 
                mmffVariant='MMFF94s'
                )
            for cid in range(frag.GetNumConformers()):
                ff = AllChem.MMFFGetMoleculeForceField(
                    frag, 
                    mp, 
                    confId=cid
                    )
                e = ff.CalcEnergy()
                res.append((cid, e))
            sorted_res = sorted(res, key=lambda x:x[1])
            rdMolAlign.AlignMolConformers(frag)
            for cid, e in sorted_res:
                frag.SetProp('CID', str(cid))
                frag.SetProp('Energy', str(e))
                w.write(frag, confId=cid)
            w.close()

        prep_log.log()

    prep_log.close()

    opt_batch_file.close()
    esp_batch_file.close()
    psi4_batch_file.close()

    del opt_batch_file
    del esp_batch_file
    del psi4_batch_file

    return conf_list