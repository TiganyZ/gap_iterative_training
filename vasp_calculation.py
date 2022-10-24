#!/usr/bin/python3
from ase import Atoms
from ase.io import read, write
from ase.calculators.vasp import Vasp
import os, shutil, subprocess
from utils import Utils
from calculations import Calculation, CalculationUtils, CalculationContainer

class VaspCalc( Calculation ):

    def __init__(self, args):
        self.name = "VaspCalc"
        self.args = args
        self.utils = Utils()
        self.calc_utils = CalculationUtils()
        self.result = {}

        self.structure = self.calc_utils.get_structure(args)


    def setup(self):
        # Copy gap files from directory to where the calculation is
        self.utils.check_keys(self.args)

        pot_dir = self.args["potential_directory"]
        out_dir = self.args["output_directory"]
        input_dir = self.args["input_directory"]

        self.utils.check_copy_file(pot_dir, "POTCAR", out_dir)

        # Will create INCAR later
        # self.utils.check_copy_file(input_dir, "INCAR", out_dir)

        shutil.copytree( f"{input_dir}/", f"{out_dir}/" )

        self.cwd = os.getcwd()
        os.chdir( out_dir )

        self.create_run_environment()


    def create_run_environment(self):
        self.utils.check_keys(self.args, keys=( "ncores", "binary" ) )

        binary = self.args["binary"]
        ncores = self.args["ncores"]
        # Make the run_vasp for the number of cores that we want
        with open("run_vasp.py", 'w') as f:
            f.write(f"""
import os
exitcode = os.system('srun -n {ncores} {binary}')
""")

        cwd = os.getcwd()

        os.environ["VASP_SCRIPT"]=f"run_vasp.py"





    def calculate(self):
        # Create the input file for the directory and then compute
        self.utils.check_keys(self.args, keys=( "structure", "input_args" ) )

        self.structure.calc = Vasp( **self.args["input_args"] )
        print(self.structure.calc)
        res =  self.structure.get_potential_energy()

        print("VASP ENERGY RESULT: ", res)
        self.result = {"result":res}#self.calc_utils.run(self.structure, self.args)

        return self.result


    def get_data(self):
        pass
