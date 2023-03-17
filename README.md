# [WIP] Adding things, fixing things

* Better support for different queuing systems
* Rdkit conformer pipeline
* Psi4 support
* Better handling of complicated charge constraints


## FragResp

FragResp is a command line tool (in theory also usable through python, please see below) that facilitates the semi-automated calculation of RESP partial charges with minimized conformation bias. The key steps are: 

* Automatic decomposition of molecules into small pieces based on a set of user-defined decomposition rules. These small pieces are called "fragments". The term is not ideal, since it is also used in other contexts of molecular sciences with varying meaning. However, we want to stick with it for now
* The fragments are capped and stored into a database
* Conformers of the fragments are generated automatically using either `openeyetoolkits` (you need a license for that) or `openbabel` (no license needed).
* Input files for Gaussian09 are generated automatically for b3lyp-level optimization and HF-level ESP calculation along with bash scripts for submission in a HPC environment.
* The convergence of the high-level calculations is checked along with other checks (bond breaking, frequencies).
* In the final step, RESP input files are generated that satisfy the fitting constraints introduced through the initial decomposition step. The program `resp` (from the AmberTools package) is invoked and partial charges are generated.

## Installation
1. Clone the repository: `git clone https://github.com/wutobias/fragresp`

2. The easiest way to use FragResp is to install it in a dedicated conda environment. Create a conda environment containing all dependencies:
`conda env create -f fragresp.yml`

3. Activate the python environment `conda activate fragresp`.

4. In the final step, install the FragResp package using either `make` or `python setup.py install` and you're done.

## Usage

### Command line

FragResp is intended to be a command line program. The executable is called `run_fragresp.py` and an overview of the expected command line arguments together with a short description of each argument can be obtained through the option `--help`.

### Python
In theory one can use FragResp from Python, but this is not really optimized. It is possible though. If you want try it, have a look at the python scripts in `fragresp/scripts/`. These scripts are invoked through the command line version of the program when a user specifies the different run modes (e.g. decomposition or resp) through the top-level routines in `fragresp/scripts/run_fragresp.py`

## Examples

Two examples can be found in `examples/`.

## Known bugs

* If a decomposition operation results in two or more identical fragments, no reasonable partial charge assignment is possible currently. For instance, consider the molecule Diethylether `CCOCC`, which is decomposed into `CCO` and `CC`. With appropriate capping groups, this results in two identical `CCOC` fragments.
* Conformer generation with openeye is not really reliable right now. This might be due to non-ideal input parameters in the conformer generation routines of FragResp.
* RDkit sometimes generates tons of warnings upon reading mol2 files, these can be suppressed.

## Things to improve

* Fix the bugs above :-)
* The defitinion of linker units in the linker input file is somewhat complicated. I think this could be largly improved by using SMIRKS patterns instead of SMILES and atom (integer) identifiers.
* Rewrite the method `resp_it_db` in `fragresp/resp_fit.py`. It is overly lengthy and contains redundant blocks of code.
* It would be nice to be able to do any step along the charge calculation workflow on any subset of molecules/fragments that is stored in the data base. Currently, one can only do it for the whole data base.
* Find a reliable way to deal with tautomers.
* Support additional QM codes such as Psi4 or Orca.
* Right now, all molecules must be in mol2 format. Support for other formats (mol, sdf, pdb, cdf) is really necessary, especially Smiles and InChI!
* Make the database searchable
* In case we run into unresolvable charge constraints, spit out a warning.
* Find a better format so store the database. Currently it is just a pickled instance of `fragresp:data_strc:fragment_db`. It would be much nicer to have a pickled dictionary as the exported version of the database.
* Currently each conformer of a given molecules contributes equally to the partial charges of each atom. This is not ideal and some smart weighting strategy should be installed. Boltzmann weights from gas phase energies are not ideal, but would be easy to implement. :-)
* Add options to allow `--capsum`, `--bycharge` or `--fittosurr` only to be applied to specific molecules/fragments. Currently it is 'all or nothing' which is not very useful for large databases.

## Final remarks

The code in this repository was written through the course of my PhD studies in the research group of Gerhard Klebe at UMR (agklebe.de). Furthermore, my former undergraduate student Janik Hedderich (now PhD student at UMR with Peter Kolb) tested the program extensivly and contributed to the test case on the FreeSolv database (see `examples/freesolv-subset`). Furthermore, Janik contributed a very nice logo :-)