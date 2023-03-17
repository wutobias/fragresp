#!/usr/bin/env python3

def main():

    import sys
    import os
    import random
    import subprocess

    NPROC = sys.argv[-2]
    COM = sys.argv[-1]

    CALCPATH = os.path.abspath(COM)
    CALCPATH = os.path.dirname(CALCPATH)
    COM      = os.path.basename(COM)

    slurm_str = f"""#!/bin/bash

#PBS -l nodes=1:ppn={NPROC}
#PBS -l walltime=24:00:00
#PBS -q home-hopper
#PBS -S /bin/bash
#PBS -N gaus.s
#PBS -j oe
#PBS -o \\${{PBS_JOBNAME}}.out
#PBS -V
        
module load gaussian
export GAUSS_SCRDIR=\\${{TMPDIR}}

cd {CALCPATH}

g09 {COM}

"""

    RANDINT = random.randint(0, 99999)
    TMPPATH = f"/tmp/g09-{RANDINT}.sh"

    with open(TMPPATH, "w") as fopen:
        fopen.write(slurm_str)

    subprocess.call(
        ["sbatch", TMPPATH]
        )

if __name__ == "__main__":

    import sys
    if len(sys.argv) != 3:
        print(
            "Usage tscc.py <NPROC> <INPUT-COM>"
            )
        exit(0)

    main()