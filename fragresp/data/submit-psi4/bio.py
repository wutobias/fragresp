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

#SBATCH -o ./psi4.out.%%j
#SBATCH -e ./psi4.err.%%j
#SBATCH -D ./
#SBATCH -J psi4.s
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task={NPROC}
#SBATCH --mem=4000
#SBATCH --partition=s.bio
#SBATCH --mail-type=none
#SBATCH --time=24:00:00

source /u/tohuefne/init_conda.sh
conda activate psi4

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}

cd {CALCPATH}

export PSI_SCRATCH=\\${{TMPDIR}}

run_psi4 {COM}

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