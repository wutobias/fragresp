## Examples

In here you find two examples that came up along my PhD studies. The first is on a series of thrombin ligands. In this series, each ligand is a combination of 2-4 amino acid fragment molecules each linked by either a amide group or sulfonamide group. Each of these two linker units is defined in the corresponding `linker.fragresp` file. Note, that also a secondary amide group is defined in order to capture proline fragment molecules.
The second example is about a subset of the FreeSolv database. This database is a nice collection of experimental small molecule solvation thermodynamics data. In addition to the currated data, you also get input files in different standard formats. For these molecules, we generated 52 linker units with the intend to cover the chemical space of the FreeSolv subset (in fact, it covers even more than that but we tested it only for this subset). **Note, that currently this second example is broken!!!** I carried out some changes to the code that turned out to be not compatible with the linker definitions in `linker.fragresp`. However, the first example on the thrombin example works perfectly fine.

### General workflow

Both examples follow the same workflow, which is demonstrated in the following. It should be transferable to both examples.
In order to generate a database with fragment molecules, three mendatory files are needed: **(A)** a linker file containing all linker units and their specificiations (explanation of what a linker actually is, follwos in a minute!); **(B)** a molecule file containing a name and file path to each molecule that should be decomposed, currently all molecules must be in mol2 format; **(C)** a input file that specifies the paths for storing data from the QM calculations and the data base. A word on how the input of these files must look like:

##### (A) The linker file

This is the most important file. It defines how molecules are becoming decomposed, how broken bonds are capped and where charges are going to be constrained during RESP charge fit. In FragResp each molecule is treated as combination of fragments and linkers. Without any restrictions towards which one (fragments or linkers) is bigger, a molecule always contains fragments and these fragments are linked to each other by linkers. That means that any molecule is something like `F1-L1-F2-L1-F3-L2-F4`, where `{F1,F2,F3,F4}` is a set of distinguishable fragments and `{L1,L2}` is a set of distinguishable linkers. This idea is well-suited for polymer molecules such as proteins or DNA, but often also works for drug molecules. In very practical sense, linkers are supposed to be ubiquitous functional groups in molecules (e.g. amide groups in proteins). But also note that linkers will have the same partial charge across all molecules, so it might not be physically reasonable for every functional group to be part of a linker.

`init/end` (string): This starts/ends a new block specifying a linker unit.

`name` (string): The name of the linker. Something descriptive is a good idea!

`terminal` (string): Indiciates whether or not this linker unit should be treated as a terminal linker. In most cases, linkers will be non-terminal. A terminal linker is one that connects only to a single fragment (e.g. a chlorine atom). One should use terminal linkers cautiously, since it runs the risk to get unrealistic partial charges for terminal atoms. This parameter can take values `on` or `off`.

`terminal_anc` (integer list): If the parameter `terminal` is set to `on`, then this parameter must be provided. It is defined similar to `lanc`/`ranc` and specifies the terminal anchor atoms, which are the atoms that are being treated as the terminal part of the linker unit.

`linker` (SMILES): This is the pattern (expressed as SMILES) that defines our actual linker. It is usually very helpful to also have relevant atoms adjacent to the actual functional group present in that SMILES pattern. For example, in the case of an amide group, the SMILES string `NC(=O)` won't help us to distinguish primary or secondary amide groups. Therefore, the actual linker must be specified as `CNC(=O)C` or as `C(C)NC(=O)C` if one wants to hit a secondary amide group.

`lcap` (SMILES): This is the SMILES representation of the left capping group. It is becoming attached to the 'left' site of the linker after a given bond connecting the left site of the linker with another fragment is broken. Note that 'left' and 'right' are in principle arbitrary, but refer to the ordering of the atoms in the linker SMILES string: the first atom in the SMILES string is on the left, the last atom on the right. This is usually equivalent to what is being capped off on the left site of the linker upon breaking bonds between linkers and fragments. However, it does not have to be equivalent. The `linker` SMILES is only for the definition of the substructure to be matched in the molecule, whereas the capping groups is what determines the fragments' actual composition. The capping groups should resemble the part that is cut off (i.e. preserving local electronic properties) but should be sufficiently small.

`lanc` (integer list): This is the left anchor of the linker and defines the atoms where the linker is 'anchored' at the fragment. The atom indices specified with this parameter correspond to the atom ordering in the SMILES string of the `linker` parameter. The anchor atoms will be treated as part of the fragment (here, the 'left' part of the fragment) and therefore all linker atoms that are not part of the anchor atoms, will be deleted and replaced by the capping group (see `lcap`).

`lancmap` (integer list): In the original implementation of the RESP procedure, the partial charges for some atoms in the linker are taken from a reference molecule (NME-ALA-ACE). This reference molecule is called 'surrogate' within FragResp and the mapping of the partial charges on those surrogate molecules on the actual linker atoms are carried out through this parameter `lancmap`. This parameter holds a list of atom indices with the same length as `lanc`. This time, the atom indices refer to the SMILES code definition of the surrogate molecule defined under `surrogate` (see below). For each atom defined in `lanc`, the partial charge will be taken from the corresponding atom in the surrogate molecule. A value of -1 in `lancmap` indicates that this atom will not be mapped from the surrogate and is fitted without contraint (except for chemical equalization). However note, that it is not strictly necessary to use a surrogate molecule in the fitting. FragResp supports partial charges without surrogate molecules, which means that the parameter `lancmap` is not needed in that case.

`rcap`, `ranc` and `rancmap` are working exactly the same as their 'left' counterpart.

`breakring (string)`: Indicates if a bond in or directly adjacent to a ring is allowed to be broken during decomposition. Parameter can be set to either `on` or `off`.

`surrogate` (SMILES): This molecule is used for deriving the charges for the cap moieties (see `lcap` and `rcap`) as well as parts of the linkers (see `lancmap` and `rancmap`) if requested.

`s_lanc` (integer list): This list of atom indices indicate the position in the surrogate molecule that are going to be replace by the left capping group specified under `lcap`.

`s_ranc` (integer list): similar to `s_lanc`

Example:
```
init
name Primary Peptide
linker CNC(=O)C
terminal off
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
```

##### (B) The molecule file

`init/end mol2` (string): Start/end of a block containing mol2 files. Each block consists of two columns, the first column holds the names of the molecules (e.g. a pdb code) and the second column holds the (relative) path to the molecule. Currently only mol2 format is supported.

Example:
```
init mol2
1K21 mol2/1K21.mol2
2ZC9 mol2/2ZC9.mol2
2ZDA mol2/2ZDA.mol2
2ZDV mol2/2ZDV.mol2
2ZF0 mol2/2ZF0.mol2
end mol2
```

##### (C) The input file

`database` (string): This is the (relative) path to the file containing the fragment database. It really is only a pickled instance of `fragresp:data_strc:fragment_db`.

`frag_dir` (string): This is the (relative) path to the directory where we will save all conformers and output from the QM calculations on the fragment molecules.

`surr_cap_dir` (string): This is the (relative) path to the directory where we will save all conformers and output from QM calcs on the capped surrogate molecules.

`mol_dir` (string): This is the (relative) path to the directory where we will save all conformers and output from QM calcs on the molecules that have not been decomposed. For these molecules we carry out a standard RESP fit.

Example:
```
database     frag_db.pickle
frag_dir     frag_dir
surr_cap_dir surr_cap_dir
mol_dir      mol_dir
```

### 1.) Generating the database

The first step is always to generate the database. In order to do so, we will execute `run_fragresp` in its `decompose` mode. Use the command:

`run_fragresp --input input.fragresp --mol2list mol2list.fragresp --linker linker.fragresp --mode decompose`

Note, that the two command line arguments `--input` and `--mode` must be present for any valid operation of `run_fragresp`. This command reads the linkers in `linker.fragresp` and applies them in order to decompose the molecules in `mol2list.fragresp`. The additional argument `--add`

### 2.) Inspect the database

In order to examine the obtained fragment molecules, perform a data dump on the whole database. Use the following command:

`run_fragresp --input input.fragresp --mode datadump`

This will create a folder `datadump`, which contains the files `frag.dat`, `mol.dat` and `surr.dat` containing information about the fragments, the non-decomposed molecules and the capped surrogate molecules, respectively. Furthermore, there is a small png under `datadump/png` for each fragment, molecule or surrogate.

### 3.) Generate conformers

After having inspected the database, we will generate conformations and input files for Gaussian09 calculations. For conformor generation we can use either `openeye` tools or `openbabel` tools, specfied by either `--pipeline openeye` or `--pipline openbabel`. Currently, the openeye pipeline is a bit more error prone than the openbabel one. Therefore, it is recommended to use openbabel. The conformer generation takes two hyperparameters `--percentage` and `--limit`. The first one, `--percentage` indicates the percentage of input conformers to select as the electrostatically least-interacting set in the first step of the [ELF](https://docs.eyesopen.com/toolkits/python/quacpactk/molchargetheory.html?) algorithm and the latter one, `--limit` indicates the number of diverse conformers to select from the least interacting conformers set in the second step of the ELF algorithm. The `--queue` argument will let FragResp write bash submissions scripts formated for different queing systems. Here we use `none`, which will not result in any preformatting and lets us run the bash script on our work station.

`run_fragresp --input input.fragresp --mode conformers --pipeline openbabel --queue none`

### 4.) Run geometry optimization

In each of the three folders `surr_cap_dir`, `mol_dir` and `frag_dir` we will find a file `submit_opt_g09.sh`. Run this bash script on your workstation or HPC resource. Make sure that Gaussian09 is installed and configured properly and that it is in your `$PATH`.

### 5.) Check optimized geometries

After the Gaussian09 runs are completed, we will perform some sanity checks on the final results of these calculations. Use the command:

`run_fragresp --input input.fragresp --mode check_opt`

This will look through the optimized structures and check if the optimized molecule is still the same as the input molecule, if the found energy minimum does not correspond to imaginary frequencies or if there is any other error. In each of the aforementioned three folders, there will be a `check_opt.log` file. Any calculation with critical output will be flagged with a `Warning !`. In case there is such a flag on a conformer of a given fragment, examine the Gaussian output file and adjust the Gaussian input parameters in the g09 input file. Some times it also helps to look at the convergence plot `opt_conv.png`, found in each of the `fragXXX` folders.

### 6.) Run ESP calculation

After we are satisfied with the optimized geometries, we want to calculate the ESP centers. This is usually less error prone and works smoothly if the geometry optimization was succesful. In order to run the ESP calculations, go to each of the three folders and run the bash script `submit_esp_g09.sh`.

### 7. Check ESP centers

After the ESP fit centers have been calculated by Gaussian09, perform some sanity checks on the results. Use the following command:

`run_fragresp --input input.fragresp --mode check_esp`

This will result in a file `check_esp.log` in each of the three folders. In case there was a problem, you will see a `Warning !` flag.

### 8. Perform the RESP fit

In the last step, we will carry out the actual RESP fitting. In the background, `run_fragresp` will always call the `resp` program  as well as other packages from the AmberTools package. So all that `run_fragresp` really does, is writing large input files for `resp` (which can be tricky) and transforming files from one format to another. In order to perform an ordinary multimolecule-multiconformer (caution, no weighting on the conformers!) resp fitting, run:

`run_fragresp --input input.fragresp --mode resp --remap --bycharge --fittosurr`

The `--remap` argument will tell `run_fragresp` to map the partial charges on the fragment molecules back to the initial set of molecules at the end. The final molecules will then be found in `remap_dir`. The `--bycharge` argument will perform a seperate resp fit for molecules with equal total charge. This is often beneficial, since it prevents the 'leaking' of charge between molecules of differing total charge (if fitted simultaneously). The `--fittosurr` argument ensures that the charges on the capping groups are taken from the surrogate molecule. Note that this is not strictly necessary and one can also fit the partial charges on the capping group during the resp fit. Another interesting argument in this context is `--capsum`, which uses slightly different boundary conditions on the charges of the capping groups. In the regular resp fit, the sum of the charges on each capping group must be equal to zero. With the `--capsum` one only requires that the sum of the partial charges of a given pair of capping groups (e.g. ACE and NME) equals to zero.
