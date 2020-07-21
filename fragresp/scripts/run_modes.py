from collections import OrderedDict

from fragresp.scripts.run_conformers import conformers
from fragresp.scripts.run_datadump import datadump
from fragresp.scripts.run_decompose import decompose_molecules
from fragresp.scripts.run_check import check_opt
from fragresp.scripts.run_check import check_esp
from fragresp.scripts.run_resp import resp

run_modes = OrderedDict()
run_modes['decompose']  = decompose_molecules
run_modes['datadump']   = datadump
run_modes['conformers'] = conformers
run_modes['check_opt']  = check_opt
run_modes['check_esp']  = check_esp
run_modes['resp']       = resp