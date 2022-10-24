#!/usr/bin/python3
from ase import Atoms
from ase.io import read, write
from ase.calculators.vasp import Vasp
import os, shutil, subprocess
from quippy.potential import Potential
from utils import Utils
from calculations import Calculation


class GapCalc( Calculation ):

    def __init__(self, args):
        self.name = "GapCalc"
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

        # self.utils.check_copy_file(input_dir, "input", out_dir)

        self.utils.check_copy_tree(input_dir, out_dir)

        if f"{input_dir}" == f"{pot_dir}":
            print("""
            ########################################################################################################
            ###---   WARNING: Input dir is the same as the potential dir for the gap calculation. Not copying.---###
            ###---            Usually this path should be a gap_files directory                               ---###
            ########################################################################################################

            """)
        if self.utils.check_key(self.args["input_args"], "quip"):
            # Copy all the xml files to the directory from the potential directory
            self.utils.copy_only_files(pot_dir, out_dir)
        else:
            if not os.path.exists(f'{out_dir}/gap_files'):
                os.mkdir(f'{out_dir}/gap_files')
            self.utils.copy_only_files(pot_dir, f'{out_dir}/gap_files')



        self.cwd = os.getcwd()
        os.chdir( out_dir )



    def calculate(self):
        # Create the input file for the directory and then compute
        self.utils.check_keys(self.args, keys=( "structure", "input_args" ) )

        if self.utils.check_key(self.args["input_args"], "quip"):
            self.setup_quip()
        else:
            self.setup_turbogap()

        self.result = self.calc_utils.run(self.structure, self.args)

        return self.result


    def setup_quip(self):
        if self.utils.check_key(self.args, "system"):
            pot_file = self.utils.check_file_dir_subdir(f'{self.args["system"]}.xml')
            gap = Potential(param_filename=pot_file)
        else:
            print(f"""
            ##############################################################################################
            ###---   WARNING: There is no system specified for this GAP calculation. Rectify this   ---###
            ##############################################################################################
            """)
            exit(1)
        # Read the configurations
        self.structure.set_calculator(gap)

    def setup_turbogap(self):
        pass

    def get_data(self):
        pass
