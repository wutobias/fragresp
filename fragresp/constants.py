from rdkit import Chem

StandardTemp    = 298.15
hartree_to_kcal = 627.50960803
kb_kcal         = 0.0019872041
kbT_kcal        = kb_kcal*StandardTemp

CHI_TETRAHEDRAL_CW  = Chem.ChiralType.CHI_TETRAHEDRAL_CW
CHI_TETRAHEDRAL_CCW = Chem.ChiralType.CHI_TETRAHEDRAL_CCW

max_confs = 1000