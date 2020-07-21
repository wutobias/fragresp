
"""
This module contains classes, methods and datastructures for 
automated derivation of RESP charges. It automates the process
of molecular decomposition, blocking group assignment, input
file generation for QM software (currently only Gaussian09) atomic
charge equalization and rebuilding of the molecule.

Written by
Tobias Hüfner
UCSD, Skaggs School of Pharmacy
and University of Marburg, Institute for Pharmaceutical Chemistry
"""
__author__     = "Tobias Hüfner"
__license__    = "MIT"
__maintainer__ = "Tobias Hüfner"
__email__      = "tobias.wulsdorf@gmail.com"

#from fragresp import decomposition
#from fragresp import data_strc
#from fragresp import gaussian_utils
#from fragresp import utils
#
#from decomposition import decompose
#from data_strc import fragment_db, connector
#from gaussian_utils import check_opt, check_esp
