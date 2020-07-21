import os

AMBERHOME=os.getenv('AMBERHOME')
if AMBERHOME == None:
    raise Warning("Environment variable AMBERHOME not set!")

ante_exe    = AMBERHOME+"/bin/antechamber"
respgen_exe = AMBERHOME+"/bin/respgen"
espgen_exe  = AMBERHOME+"/bin/espgen"
resp_exe    = AMBERHOME+"/bin/resp"