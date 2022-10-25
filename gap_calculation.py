#!/usr/bin/python3
from ase import Atoms
from ase.io import read, write
from ase.calculators.vasp import Vasp
import os, shutil, subprocess
from quippy.potential import Potential
from utils import Utils
from calculations import Calculation, CalculationUtils, CalculationContainer


class GapCalc( Calculation ):

    def __init__(self, args):
        self.name = "GapCalc"
        self.args = args
        self.utils = Utils()
        self.calc_utils = CalculationUtils()
        self.result = {}
        self.structure = self.calc_utils.get_structure(args)


    def find_potential_file(self, directory):
        # First check this current directory

        pot_file = self.utils.check_file_dir_subdir(f'{self.args["system"]}.xml', dir='.', subdir="gap_files")

        if pot_file is None:
            pot_file = self.utils.check_file_dir_subdir(f'{self.args["system"]}.xml', dir=directory, subdir="gap_files")

        if pot_file is None:
            print("""
            ###################################################################################################################################################
            ###---   FATAL: Could not find the potential file {self.args['system']}.xml in ./ or ./gap_files or {directory} or {directory}/gap_files   ---###
            ###################################################################################################################################################

            """)
            exit(1)
        else:
            return os.path.abspath(pot_file)


    def setup(self):
        # Copy gap files from directory to where the calculation is

        out_dir = self.args["output_directory"]
        pot_file = self.find_potential_file(out_dir)
        self.pot_path = os.path.abspath(pot_file)
        if self.utils.check_key(self.args, "input_args"):
            if self.utils.check_key(self.args["input_args"], "quip"):
                self.setup_quip(pot_file, out_dir)
            else:
                self.setup_turbogap()
        else:
            self.setup_turbogap()



    def run(self):
        # Create the input file for the directory and then compute
        self.result = self.calc_utils.run(self.structure, self.args)

        return self.result


    def setup_quip(self, pot_file, out_dir):
        print(">>>   Setting up QUIP gap calculation <<<")
        gap = Potential(param_filename=pot_file, directory = out_dir)
        self.structure.set_calculator(gap)
        print(">>>   Assigning calculator <<<")
        self.calc_func = Potential
        self.calc_args = {"param_filename":self.pot_path, "directory" : out_dir}



    def setup_turbogap(self):
        pass

    def get_data(self):
        pass
