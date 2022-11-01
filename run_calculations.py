#!/usr/bin/python3

from calculations import CalculationData, CalculationContainer
from neb_interface import NEB_data, NEB_interface
import os, shutil, subprocess, pickle
from ase.io import read, write


# These are actually the results of a relaxed gap calculation
files = [f"C+Br_GapCalc_Br_graphene_before_relaxation.xyz", f"C+Br_GapCalc_Br_graphene_before_relaxation_translated.xyz"]


base = files[0].split('.')[0]
istructure = read(files[0])
fstructure = read(files[1])


gap_input_args = {"quip" : True,
                  "pot1" : "carbon.xml",
                  "pot2" : "CBr.xml"}


args = CalculationData( binary              = "-",
                        input_directory     = "./",
                        potential_directory = "./gap_files",
                        output_directory    = f"output_NEB_{base}",
                        structure           = istructure,
                        input_args          = gap_input_args,
                        ncores              = 1,
                        system              = "C+Br", 
                        run_calc_type       = "energy"
                       )

# Do neb calculation using gap


from gap_calculation import GapCalc

calculation_method = GapCalc

c1 = CalculationContainer(calculation_method,  args )



args = CalculationData( binary              = "-",
                        input_directory     = "./",
                        potential_directory = "./gap_files",
                        output_directory    = f"output_NEB_{base}",
                        structure           = fstructure,
                        input_args          = gap_input_args,
                        ncores              = 1,
                        system              = "C+Br", 
                        run_calc_type       = "energy"
                       )

c2 = CalculationContainer(calculation_method,  args )

neb_args = NEB_data(images = [c1.method, c2.method],
                    n_images = 5,
                    climb = False
                    )

with open("neb_args.pkl", "wb") as f:
    pickle.dump( neb_args, f)

n = CalculationContainer(NEB_interface, neb_args )
n.method.name = f"{n.method.name}_{base}"
n.run()



# print(f"\n -> Directory <-\n  {os.getcwd()} \n\n")
# # Now create DFT calculation with this structure

# # Final trajectory is in the path

# shutil.copy(f"{c.method.path}/{args.run_calc_type}_{c.method.name}.xyz", f"current_calc_GAP_{base}.xyz")

# # Now pass on to vasp

# # Now have completed vasp calculaton
# print(f"file = \"current_calc_GAP_{base}.xyz\"")

# from utils import Utils
# u = Utils()

# vfile = "run_vasp_of_gap.py"
# rstr=f"current_calc_GAP_{base}.xyz"

# output = "vasp_calculations_to_run.sh"
# command = f"python3.9 {vfile} {rstr}\n"

# if os.path.exists(output):
#     mode = 'a'
# else:
#     mode = 'w'

# with open(output, mode) as f:
#     f.write( command )

