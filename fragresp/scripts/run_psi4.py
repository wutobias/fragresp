#!/usr/bin/env python3

def read_cominput(inputcom):

    import numpy as np

    nproc = 8
    mem   = 2
    input_str = ""

    with open(inputcom, "r") as fopen:
        found_input_line = False
        reading_comment = 0
        for line in fopen:
            line = line.rstrip().lstrip()
            if reading_comment == 0:
                if line.startswith("!"):
                    continue
                if len(line) == 0:
                    continue
                if line == "\n":
                    continue

                segments = line.split()
                if segments[0].startswith("%%nproc"):
                    nproc = line.split("=")[-1]
                if segments[0].startswith("%%mem"):
                    mem = line.split("=")[-1]
                    if mem.endswith("MB"):
                        mem = mem.replace("MB", "")
                        mem = float(mem)
                        mem = np.ceil(mem/1000.)
                    elif mem.endswith("GB"):
                        mem = mem.replace("GB", "")
                        mem = float(mem)
                    else:
                        raise ValueError(
                            "mem keyword not understood"
                            )
                if line.startswith("#"):
                    reading_comment = 1

            elif reading_comment == 1:
                reading_comment = 2
            elif reading_comment == 2:
                reading_comment = 3
            elif reading_comment == 3:
                reading_comment = 4
            elif reading_comment == 4:    
                if line.startswith("!"):
                    break
                if len(line) == 0:
                    break
                if line == "\n":
                    break
                input_str += line + "\n"

    return input_str, nproc, mem


def run_esp(input_str):

    import qcelemental as qcel
    from openff.toolkit.topology import Molecule
    from pint import UnitRegistry
    import os
    from openff.recharge.esp import ESPSettings
    from openff.recharge.esp.psi4 import Psi4ESPGenerator
    from openff.recharge.grids import MSKGridSettings

    ureg = UnitRegistry()

    qc_data_settings = ESPSettings(
        method="b3lyp", 
        basis="6-31G*", 
        grid_settings=MSKGridSettings()
    )

    qcemol   = qcel.models.Molecule.from_data(input_str)
    geometry = qcemol.geometry * ureg.bohr

    offmol = Molecule()
    offmol.add_atom(
        int(qcemol.atomic_numbers[0]),
        int(qcemol.molecular_charge),
        False
        )
    nat = len(qcemol.atomic_numbers)
    for idx in range(1, nat):
        offmol.add_atom(
            int(qcemol.atomic_numbers[idx]),
            0,
            False
            )

    if "PSI_SCRATCH" in os.environ:
        conformer, grid, esp, electric_field = Psi4ESPGenerator.generate(
            offmol, 
            geometry.to(ureg.angstrom), 
            qc_data_settings, 
            minimize=True,
            directory=os.environ["PSI_SCRATCH"]
        )
    else:
        conformer, grid, esp, electric_field = Psi4ESPGenerator.generate(
            offmol, 
            geometry.to(ureg.angstrom), 
            qc_data_settings, 
            minimize=True,
        )

    return grid, esp

def main():

    import sys
    if len(sys.argv) != 2:
        print(
            "Usage run_qm.py <input.com>"
            )
        exit(0)

    inputcom = sys.argv[-1]

    input_str, nproc, mem = read_cominput(inputcom)
    grid, esp = run_esp(input_str)


def entry_point():

    main()

if __name__ == '__main__':

    entry_point()