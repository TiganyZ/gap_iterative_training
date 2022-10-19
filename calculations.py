#!/usr/bin/python3

from abc import ABC, abstractmethod
from typing import Type, TypeVar
from datetime import datetime
from ase import Atoms
from ase.io import read, write
from ase.calculators.vasp import Vasp
import os, shutil, subprocess
import quippy
from quippy.potential import Potential

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
        os.environ["VASP_SCRIPT"]=f"{cwd}/run_vasp.py"


    def get_energy(self):
        # Create the input file for the directory and then compute
        self.utils.check_keys(self.args, keys=( "structure", "input_args" ) )

        structure = self.args["structure"]
        structure.calc = Vasp( **self.args["input_args"] )
        energy = structure.get_potential_energy()

        return energy


    def get_data(self):
        pass

class GapCalc( Calculation ):

    def __init__(self, args):
        self.name = "GapCalc"
        self.args = args
        self.utils = Utils()

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


    def get_energy(self):
        if self.utils.check_key(self.args["input_args"], "quip"):
            energy = self.get_energy_quip()
        else:
            # Use turbogap
            energy = 0.0

        return energy

    def get_energy_quip(self):
        if self.utils.check_key(self.args, "system"):
            gap = Potential(param_filename=f'{self.args["system"]}.xml')
        else:
            print(f"""
            ##############################################################################################
            ###---   WARNING: There is no system specified for this GAP calculation. Rectify this   ---###
            ##############################################################################################
            """)
            exit(1)
        # Read the configurations
        structure = self.args["structure"]
        structure.set_calculator(gap)
        energy = structure.get_potential_energy()

        print(f">>> GAP with QUIP: energy  = {energy} <<< ")
        return energy


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

    test = "Vasp"
    test = "Gap"


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
                "structure"           : structure,
                "input_args"          : gap_input_args,
                "ncores"              : ncores,
                "system"              : system

        }

        calculation_method = GapCalc

        c = CalculationContainer(calculation_method, args )

        c.run()


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


        args ={ "binary"              : binary,
                "potential_directory" : potential_directory,
                "input_directory"     : input_directory,
                "output_directory"    : output_directory,
                "structure"           : structure,
                "input_args"          : vasp_input_args,
                "ncores"              : ncores

        }

        calculation_method = VaspCalc

        c = CalculationContainer(calculation_method, args )

        c.run()
