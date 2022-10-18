#!/usr/bin/python3

from abc import ABC, abstractmethod
from typing import Type, TypeVar
from datetime import datetime
from ase import Atoms
from ase.io import read, write
from ase.calculators.vasp import Vasp
import os, shutil, subprocess

from utils import Utils
# from process_calculation import GAP_to_VASP, VASP_to_GAP


class Calculation(ABC):

    @abstractmethod
    def setup():
        pass

    @abstractmethod
    def get_energy():
        pass

    @abstractmethod
    def get_data():
        pass



class VaspCalc( Calculation ):

    def __init__(self, args):
        self.name = "VaspCalc"
        self.args = args
        self.utils = Utils()


    def setup(self):
        # Copy gap files from directory to where the calculation is
        self.utils.check_keys()

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


    def create_run_environment():
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
        os.environ["VASP_SCRIPT"]=f"{cwd}/run_vasp.py"


    def get_energy():
        # Create the input file for the directory and then compute
        self.utils.check_keys(self.args, keys=( "structure", "input_args" ) )

        structure = self.args["structure"]
        structure.calc = Vasp( **self.args["input_args"] )
        energy = structure.get_potential_energy()

        return energy


    def get_data():

        return


class GapCalc( Calculation ):

    def __init__(self, args):
        self.name = "GapCalc"
        self.args = args

    def setup(self):
        # Copy gap files from directory to where the calculation is
        self.utils.check_keys()

        pot_dir = self.args["potential_directory"]
        out_dir = self.args["output_directory"]
        input_dir = self.args["input_directory"]

        # self.utils.check_copy_file(input_dir, "input", out_dir)

        shutil.copytree( f"{pot_dir}", f"{out_dir}/" )
        shutil.copytree( f"{input_dir}/", f"{out_dir}/" )

        self.cwd = os.getcwd()
        os.chdir( out_dir )


    def get_energy():
        self.utils.check_keys(("system",))
        gap = Potential(param_filename=f'{self.args["system"]}.xml')

        # Read the configurations
        structure = self.args["structure"]
        structure.set_calculator(gap)
        energy = structure.get_potential_energy()

        return energy


    def get_data():

        return





class CalculationContainer:
    def __init__(self, calculation_method: Type[Calculation],
                 binary: str, potential_directory: str, output_directory: str,
                 structure: Atoms, input_args: dict):


        self.binary = binary
        self.potential_directory = potential_directory
        # self.images = sorted(os.listdir(self.potential_directory))

        args = { "binary": binary,
                 "input_directory": input_directory,
                 "potential_directory": potential_directory,
                 "structure": structure,
                 "input_args": input_args
                }

        self.method = calculation_method(args)
        self.method_name = self.method.name

        now = datetime.now()
        self.dt = now.strftime("%Y-%m-%d--%H-%M-%S")

        self.output_directory = self.create_output_directory()
        self.method.args["output_directory"] = self.output_directory

        self.run_setup = False
        self.run_energy = False
        self.run_get_data = False


    def create_output_directory(self):
        dir_name = f"{self.method_name}_{self.dt}"
        if not os.path.exists(dir_name):
            os.mkdir(dir_name)
        return dir_name


    def run(self):

        self.method.setup()
        self.run_setup = True

        self.method.get_energy()
        self.run_energy = True

        self.method.get_data()
        self.run_get_data = True


if __name__ == "__main__":


    input_directory = "./input_dir"
    output_directory = "./output_dir"
    potential_directory = "./gap_files"

    binary = "turbogap"

    structure = read("POSCAR", format="vasp")

    vasp_input_args = {"prec" : 'Accurate',
                        "xc" : 'PBE',
                        "ibrion" : = -1,
                        "algo" : = 'Normal',
                        "potim" : = 1.0,
                        "ediff" : = 1e-6,
                        "lwave" : False,
                        "lcharg" : False,
                        "ncore" : 16,
                        "kpar" : 4
    }


    calculation_method = GapCalc

    c = CalculationContainer(calculation_method  = calculation_method,
                             binary              = binary,
                             potential_directory = potential_directory,
                             output_directory    = output_directory,
                             structure           = structure,
                             input_args          = vasp_input_args
                             )

    c.run()
