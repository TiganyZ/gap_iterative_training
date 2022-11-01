#!/usr/bin/python3
from abc import ABC, abstractmethod
from typing import Type, TypeVar, Union
from datetime import datetime
from ase import Atoms
from ase.io import read, write
from ase.calculators.vasp import Vasp
import os, shutil, subprocess
from utils import Utils
from ase.optimize import BFGS
from json import JSONEncoder

# from process_calculation import GAP_to_VASP, VASP_to_GAP

from dataclasses import dataclass, field

@dataclass
class CalculationData:
    input_directory: str
    output_directory: str
    potential_directory: str
    binary: str
    input_args: dict
    structure: Atoms
    ncores: int
    system: str
    make_dirs: bool = True
    driver_args: dict = field(default_factory=dict)
    batch: bool = False
    run_calc_type: str = "energy"
    run_calc_args: Union[dict, None]= None
    result: Union[dict, None]= None


class CalculationDataEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, Atoms):
            return obj.__dict__
        # Base class default() raises TypeError:
        return JSONEncoder.default(self, obj)





class Calculation(ABC):

    @abstractmethod
    def setup():
        pass

    @abstractmethod
    def run():
        pass

    @abstractmethod
    def save():
        pass

    @abstractmethod
    def save_state():
        pass

    @abstractmethod
    def get_data():
        pass


class CalculationUtils:
    def __init__(self):
        self.utils = Utils()


    def determine_calculation_type(self, args):
        calc_types = [ "energy", "optimize" ]

        calc_type = args.run_calc_type
        calc_args = args.run_calc_args

        if not (calc_type in calc_types):
            print(f"""
            #################################################################################
            ###---   WARNING: {calc_type} not in the implemented calc types: Exiting   ---###
            #################################################################################
            """)
            raise ValueError
        else:
            return calc_type, calc_args


    def run(self, structure, args, name = "", path="."):

        calc_type, calc_args = self.determine_calculation_type(args)

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
            res =  self.run_calc(structure.get_potential_energy, calc_args)

        if calc_type == "optimize":
            res = BFGS( structure, trajectory = f"{path}/{calc_type}_{name}.traj" )
            self.run_calc(res.run, calc_args)
            write(f"{path}/{calc_type}_{name}.xyz", read(f"{path}/{calc_type}_{name}.traj", index=':'), format="extxyz")
            result["optimized_structure"] = read(f"{path}/{calc_type}_{name}.traj")

        result["result"] = res

        print(result)
        return result
        

    def run_calc(self, func, args):
        if args is not None:
            return func(**args)
        else:
            return func()




class CalculationContainer:
    def __init__(self, calculation_method: Type[Calculation], args):

        self.args = args
        self.method = calculation_method(self.args)
        self.method_name = self.method.name

        now = datetime.now()
        self.dt = now.strftime("%Y-%m-%d--%H-%M-%S")

        if args.output_directory == "":
            output_directory = self.create_output_directory()
            self.method.args.output_directory = output_directory
        out_dir = args.output_directory
        input_dir = args.input_directory

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
        self.method.save_state()
        self.run_setup = True

        self.method.utils.wrap_function(f"{self.method_name}.run", self.method.run, "calculating")
        self.method.save_state()
        self.run_calculation = True

        self.method.utils.wrap_function(f"{self.method_name}.save", self.method.save, "saving")
        self.method.save_state()
        self.saved_calculation = True

        self.method.args.result["method"] = self.method_name

        self.method.get_data()
        self.method.save_state()
        self.run_get_data = True


if __name__ == "__main__":

    test = "Vasp"
    test = "Gap"


    if test == "Gap":

        input_directory = "./"
        output_directory = "output_relaxation"
        potential_directory = "./gap_files"

        binary = "turbogap"
        ncores = 128

        structure = read(f"Br_graphene_before_relaxation.xyz", format="extxyz")
        #        structure = read(f"graphene_before_relaxation.xyz", format="extxyz")


        gap_input_args = {"quip" : True,
                          "pot1" : "carbon.xml",
                          "pot2" : "CBr.xml"}

        system = "CBr"
        args = CalculationData( binary              = binary,
                                potential_directory = potential_directory,
                                input_directory     = input_directory,
                                output_directory    = output_directory,
                                structure           = structure,
                                input_args          = gap_input_args,
                                ncores              = ncores,
                                system              = system, 
                                run_calc_type       = "optimize",
                                run_calc_args       = { "fmax" : 0.001 }
                               )

        from gap_calculation import GapCalc

        calculation_method = GapCalc

        c = CalculationContainer(calculation_method,  args )
        c.method.name = "CBr_Br_ads_graphene_initial"
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


        args = CalculationData( binary              = binary,
                                potential_directory = potential_directory,
                                input_directory     = input_directory,
                                output_directory    = output_directory,
                                structure           = structure,
                                input_args          = vasp_input_args,
                                ncores              = ncores
                               )

        from vasp_calculation import VaspCalc
        calculation_method = VaspCalc

        c = CalculationContainer(calculation_method,  args )

        c.run()

        print(c.method.result)
