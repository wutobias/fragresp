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

#SBATCH -o ./g09.out.%%j
#SBATCH -e ./g09.err.%%j
#SBATCH -D ./
#SBATCH -J gaus.s
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task={NPROC}
#SBATCH --mem=4000
#SBATCH --partition=s.bio
#SBATCH --mail-type=none
#SBATCH --time=24:00:00

module load gaussian/g09 
source \\${{g09root}}/g09/bsd/g09.profile

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
            "Usage bio.py <NPROC> <INPUT-COM>"
            )
        exit(0)

    main()