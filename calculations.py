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
    def run():
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
        out_dir = args["output_directory"]
        input_dir = args["input_directory"]

        self.utils = Utils()
        self.utils.check_copy_tree(input_dir, out_dir)

        self.path = os.path.abspath(out_dir)

        self.run_setup = False
        self.run_calculation = False
        self.run_get_data = False



    def create_output_directory(self):
        dir_name = f"{self.method_name}_{self.dt}"
        if not os.path.exists(dir_name):
            os.mkdir(dir_name)
        return dir_name


    def run(self):

        self.method.utils.wrap_function(f"{self.method_name}.setup", self.method.setup, "setup")
        self.run_setup = True

        self.method.utils.wrap_function(f"{self.method_name}.run", self.method.run, "calculating")
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

        from gap_calculation import GapCalc

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
                           #                            "txt": "OUTCAR_test"
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
                "ncores"              : ncores

        }

        from vasp_calculation import VaspCalc
        calculation_method = VaspCalc

        c = CalculationContainer(calculation_method,  args )

        c.run()

        print(c.method.result)
