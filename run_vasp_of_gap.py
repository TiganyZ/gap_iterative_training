#!/usr/bin/python3

from calculations import CalculationData, CalculationContainer
import os, shutil, sys
from ase.io import read, write

print(sys.argv)

file = sys.argv[1]

base = file.split('.')[0]
gap_rel_struc = read(file, format="extxyz")

# Now compute single shot DFT calculation 

vasp_input_args = {"prec" : "Accurate",
                   "ibrion" : -1,
                   "encut" : 650,
                   "ediff" : 1.0e-08,
                   "ediffg" : -1e-03,
                   "ismear" : 0,
                   "sigma" : 0.1,
                   "algo" : "fast",
                   "lwave" : False,
                   "lcharg" : False,
                   "istart" : 0,
                   "nsw" : 0,
                   "ncore" : 16 ,
                   "isif" : 2,
                   "kspacing" : 0.25,
                   "kgamma" : True,

                   'xc' : "PBE"

}


#
os.environ["VASP_PP_PATH"] = "/projappl/project_2006384/vasp/potentials"

from utils import Utils
u = Utils()

# u.piped_subprocess("echo \"$VASP_PP_PATH\"")
# u.piped_subprocess("ls -d \"$VASP_PP_PATH\"")

#    os.symlink( os.path.join(os.environ["VASP_PP_PATH"], "potpaw_LDA"), os.path.join(os.environ["VASP_PP_PATH"], "potpaw") )
# os.symlink( os.path.join(os.environ["VASP_PP_PATH"], "potpaw_PBE"), os.path.join(os.environ["VASP_PP_PATH"], "potpaw_GGA") )    

args = CalculationData( binary              = "/appl/soft/phys/vasp/6.3.0/gcc-11.2.0/bin/vasp_std",
                        potential_directory = "",
                        input_directory     = ".",
                        output_directory    = f"output_relaxation_VASP_{base}",
                        structure           = gap_rel_struc,
                        input_args          = vasp_input_args,
                        ncores              = 128,
                        system              = "CBr", 
                        run_calc_type       = "energy"
                       )

from vasp_calculation import VaspCalc
calculation_method = VaspCalc

d = CalculationContainer(calculation_method,  args )
d.method.name = f"{d.method.name}_{base}"
d.run()

print(d.method.args.result)
