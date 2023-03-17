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

    slurm_str = f"""#!/bin/bash -l

#PBS -l nodes=1:ppn={NPROC}
#PBS -l walltime=24:00:00
#PBS -q home-hopper
#PBS -S /bin/bash
#PBS -N psi4.s
#PBS -j oe
#PBS -o \\${{PBS_JOBNAME}}.out
#PBS -V

source /home/thuefner/init_conda.sh
conda activate psi4

export OMP_NUM_THREADS={NPROC}
export MKL_NUM_THREADS={NPROC}

cd {CALCPATH}

run_psi4 {COM}

"""

    RANDINT = random.randint(0, 99999)
    TMPPATH = f"/tmp/g09-{RANDINT}.sh"

    with open(TMPPATH, "w") as fopen:
        fopen.write(slurm_str)

    subprocess.call(
        ["qsub", TMPPATH]
        )

if __name__ == "__main__":

    import sys
    if len(sys.argv) != 3:
        print(
            "Usage bio.py <NPROC> <INPUT-COM>"
            )
        exit(0)

    main()