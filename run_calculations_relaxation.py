#!/usr/bin/python3

from calculations import CalculationData, CalculationContainer
import os, shutil, subprocess, sys
from ase.io import read, write


#files = [f"Br_graphene_before_relaxation.xyz", f"Br_graphene_before_relaxation_translated.xyz"]


if len(sys.argv) >= 2:
    file=sys.argv[1]
    print(f"###---   Performing GAP relaxation of {file}   ---###")
if len(sys.argv) >= 3:
    prefix=sys.argv[2]
else:
    prefix="output_relaxation_GAP"
print(f"         > Directory prefix = {prefix}")


base = file.split('.')[0]
structure = read(file, format="extxyz")


gap_input_args = {"quip" : True,
                  "pot1" : "carbon.xml",
                  "pot2" : "CBr.xml"}


args = CalculationData( binary              = "-",
                        input_directory     = "./",
                        potential_directory = "./gap_files",
                        output_directory    = f"{prefix}_{base}",
                        structure           = structure,
                        input_args          = gap_input_args,
                        ncores              = 1,
                        system              = "C+Br", 
                        run_calc_type       = "optimize",
                        run_calc_args       = { "fmax" : 0.001 }
                       )

from gap_calculation import GapCalc
calculation_method = GapCalc

c = CalculationContainer(calculation_method,  args )
c.method.name = f"{c.method.name}_{base}"
c.run()

print(c.method.args.result)

print(f"\n -> Directory <-\n  {os.getcwd()} \n\n")
# Now create DFT calculation with this structure

# Final trajectory is in the path

shutil.copy(f"{c.method.path}/{args.run_calc_type}_{c.method.name}.xyz", f"current_calc_GAP_{base}.xyz")

# Now pass on to vasp, if there are two potentials, then make two vasp calculations
# Delete the br species
element = 'Br'
atoms = read(f"current_calc_GAP_{base}.xyz")
del atoms[[atom.index for atom in atoms if atom.symbol == element]]
write(f"current_calc_GAP_{base}_no_{element}.xyz", atoms)




# Now have completed vasp calculaton
print(f"file = \"current_calc_GAP_{base}.xyz\"")

from utils import Utils
u = Utils()

vfile = "run_vasp_of_gap.py"
for rstr in [f"current_calc_GAP_{base}.xyz" , f"current_calc_GAP_{base}_no_{element}.xyz" ]:
    output = "vasp_calculations_to_run.sh"
    command = f"python3.9 {vfile} {rstr}\n"

    if os.path.exists(output):
        mode = 'a'
    else:
        mode = 'w'

    with open(output, mode) as f:
        f.write( command )

