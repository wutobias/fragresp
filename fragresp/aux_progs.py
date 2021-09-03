

### see https://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

ante_exe    = which("antechamber")
respgen_exe = which("respgen")
espgen_exe  = which("espgen")
resp_exe    = which("resp")

if ante_exe == None:
    raise Warning("antechamber not found")
if respgen_exe == None:
    raise Warning("respgen not found")
if espgen_exe == None:
    raise Warning("espgen not found")
if resp_exe == None:
    raise Warning("resp not found")