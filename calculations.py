#!/usr/bin/python3

from abc import ABC, abstractmethod
from typing import Type, TypeVar
from datetime import datetime
from ase import Atoms
from ase.io import read, write
from ase.calculators.vasp import Vasp
import os, shutil, subprocess
# from quippy.potential import Potential as QuippyPotential
import quippy

from utils import Utils
# from process_calculation import GAP_to_VASP, VASP_to_GAP


class Calculation(ABC):

    @abstractmethod
    def setup():
        pass

    @abstractmethod
    def calculate():
        pass

    @abstractmethod
    def get_data():
        pass


class CalculationUtils:
    def __init__(self):
        self.utils = Utils()


    def get_structure(self, args):

        if self.utils.check_key(args, "structure" ):
            return args["structure"]
        else:
            print("""
            ########################################################################################
            ###---   WARNING: Structure not passed to arguments. Specify before calculation   ---###
            ########################################################################################
            """)
            return None


    def determine_calculation_type(self, args):
        calc_types = [ "energy" ]

        calc_type = "energy"
        for c in calc_types:
            if self.utils.check_key(args, c):
                calc_type = c

        return calc_type


    def run(self, structure, args):

        calc_type = self.determine_calculation_type(args)

        if structure is None:
            print(f"""
            ###############################################################################
            ###---   FATAL: Structure not provided for the {calc_type} calculation   ---###
            ###############################################################################
            """)
            exit(1)

        result = {"calc_type" : calc_type}

        if calc_type == "energy":
            print(structure)
            res =  structure.get_potential_energy()


        result["result"] = res

        print(result)
        return result
        

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
        self.utils.check_keys(self.args, keys=( "ncores", "binary", "driver_args" ) )
        self.args["driver_args"]["ncores"] = self.args["binary"]
        self.args["driver_args"]["binary"] = self.args["ncores"]

        driver_template = self.utils.get_driver_template(self.args["driver_args"])

        with open("vasp_driver.sh", 'w') as f:
            f.write(driver_template)

        os.chmod("vasp_driver.sh", 0o777)





    def calculate(self):
        # Create the input file for the directory and then compute
        self.utils.check_keys(self.args, keys=( "structure", "input_args" ) )

        self.structure.calc = Vasp( **self.args["input_args"] )

        # Write the new INCAR file



        print(self.structure.calc)
        res =  self.structure.get_potential_energy()

        print("VASP ENERGY RESULT: ", res)
        self.result = {"result":res}#self.calc_utils.run(self.structure, self.args)

        return self.result


    def get_data(self):
        pass


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




class CalculationContainer:
    def __init__(self, calculation_method: Type[Calculation], args):

        self.args = args

        self.method = calculation_method(self.args)
        self.method_name = self.method.name

        now = datetime.now()
        self.dt = now.strftime("%Y-%m-%d--%H-%M-%S")

        if not "output_directory"  in args.keys():
            output_directory = self.create_output_directory()
            self.method.args["output_directory"] = output_directory


        self.run_setup = False
        self.run_calculation = False
        self.run_get_data = False


    def create_output_directory(self):
        dir_name = f"{self.method_name}_{self.dt}"
        if not os.path.exists(dir_name):
            os.mkdir(dir_name)
        return dir_name


    def run(self):

        self.method.setup()
        self.run_setup = True

        self.method.calculate()
        self.run_calculation = True
        self.method.result["method"] = self.method_name

        self.method.get_data()
        self.run_get_data = True


if __name__ == "__main__":

    test = "Vasp"
    #test = "Gap"


    if test == "Gap":

        input_directory = "input_dir"
        output_directory = "output_dir"
        potential_directory = "input_dir/gap_files"

        binary = "turbogap"
        ncores = 128

        structure = read(f"{input_directory}/POSCAR", format="vasp")

        gap_input_args = {"quip" : True}

        system = "CBr"

        args ={ "binary"              : binary,
                "potential_directory" : potential_directory,
                "input_directory"     : input_directory,
                "output_directory"    : output_directory,
                "input_args"          : gap_input_args,
                "structure"           : structure,
                "ncores"              : ncores,
                "system"              : system

        }

        calculation_method = GapCalc

        c = CalculationContainer(calculation_method,  args )

        c.run()

        print(c.method.result)


    if test == "Vasp":

        input_directory = "input_dir"
        output_directory = "output_dir"
        potential_directory = "input_dir"

        binary = "/appl/soft/phys/vasp/6.3.0/gcc-11.2.0/bin/vasp_std"
        ncores = 128

        structure = read(f"{input_directory}/POSCAR", format="vasp")

        vasp_input_args = {"prec" : 'Accurate',
                            "xc" : 'PBE',
                            "ibrion" : -1,
                            "algo" :  'Normal',
                            "potim" :  1.0,
                            "ediff" :  1e-6,
                            "lwave" : False,
                            "lcharg" : False,
                            "ncore" : 16,
                            "kpar" : 4
        }

        driver_args =  {"job-name": "vasp_quip",
                        "account": "project_2006384",
                        "queue" : "medium",
                        "ntasks-per-node": 128,
                        "nodes": 1,
                        "walltime": "0-00:10:00",
                        "output": "out_vasp_quip",
                        "modules": ["vasp/6.3.0"],
                        "export_paths": ["VASP_PP_PATH=/projappl/project_2006384/vasp/potentials"],
                        "command": "srun"}




        args ={ "binary"              : binary,
                "potential_directory" : potential_directory,
                "input_directory"     : input_directory,
                "output_directory"    : output_directory,
                "structure"           : structure,
                "input_args"          : vasp_input_args,
                "ncores"              : ncores,
                "driver_args"         : driver_args

        }

        calculation_method = VaspCalc

        c = CalculationContainer(calculation_method,  args )

        c.run()

        print(c.method.result)
