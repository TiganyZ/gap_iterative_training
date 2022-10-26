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

        if self.utils.check_key(args, "make_dirs"):
            self.make_dirs = self.args["make_dirs"]
        else:
            self.make_dirs = True


    def setup(self):
        # Copy gap files from directory to where the calculation is

        self.utils.check_keys(self.args)
        out_dir = self.args["output_directory"]
        pot_dir = self.args["potential_directory"]
        self.pot_path = os.path.abspath(f"{pot_dir}/POTCAR" )


        if self.make_dirs:
            self.copy_potential(out_dir)
            self.path = os.path.abspath(out_dir)
        else:
            self.path = os.path.abspath("./")

        self.utils.check_keys(self.args, keys=( "structure", "input_args" ) )

        self.args["input_args"]["directory"] = self.path
        self.structure.calc = Vasp( **self.args["input_args"] )
        self.calc_func = Vasp
        self.calc_args = self.args["input_args"]

        self.create_run_environment(out_dir)


    def copy_potential(self, dir):
        self.utils.check_keys(self.args)
        out_dir = self.args["output_directory"]
        pot_dir = self.args["potential_directory"]
        self.utils.check_copy_file(f"{pot_dir}/POTCAR", dir)

    def create_run_environment(self, out_dir):
        self.utils.check_keys(self.args, keys=( "ncores", "binary" ) )
        self.utils.check_keys(self.args, keys=( "ncores", "binary", "driver_args" ) )
        self.args["driver_args"]["ncores"] = self.args["binary"]
        self.args["driver_args"]["binary"] = self.args["ncores"]

        if self.utils.check_key(self.args, "batch"):
            driver_template = self.utils.get_driver_template(self.args["driver_args"])

        else:
            binary = self.args["binary"]
            ncores = self.args["ncores"]
            # Make the run_vasp for the number of cores that we want
            with open(f"{out_dir}/run_vasp.py", 'w') as f:
                f.write(f"""
import os
exitcode = os.system('srun -n {ncores} {binary}')
    """)

            cwd = os.getcwd()

            os.environ["VASP_SCRIPT"]=os.path.abspath(f"{out_dir}/run_vasp.py")





    def run(self):


        # Create the input file for the directory and then compute
        # self.utils.check_keys(self.args, keys=( "structure", "input_args" ) )

        # self.structure.calc = Vasp( **self.args["input_args"] )
        # self.calc = Vasp( **self.args["input_args"] )

        # print(self.structure.calc)
        self.result = self.calc_utils.run(self.structure, self.args)

        # print("VASP ENERGY RESULT: ", self.result)
        # self.result = {"result":self.result}#self.calc_utils.run(self.structure, self.args)

        return self.result



    def save(self, prefix=""):
        if len(prefix) == 0:
            prefix = self.name
        name = self.utils.get_save_name(self.path, self.result, prefix)
        self.structure.calc.write_json(name)



    def get_data(self):
        pass
