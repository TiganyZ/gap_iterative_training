#!/usr/bin/python3
from ase import Atoms
from ase.io import read, write
from ase.calculators.vasp import Vasp
from ase.neb import NEB
from ase.optimize import BFGS
from ase.constraints import FixAtoms
import os, shutil, subprocess, copy
from utils import Utils
from calculations import Calculation, CalculationUtils, CalculationContainer

class NEB_interface(Calculation):
    def __init__(self, args):

        self.name = "NebCalc"
        self.args = args
        self.images = args["images"]
        self.utils = Utils()

        if self.utils.check_key(args, "climb"):
            self.climb = self.args["climb"]
        else:
            self.climb = False


    def get_calc(self):
        if self.utils.check_key(self.args, "calc_func") and self.utils.check_key(self.args, "calc_args"):
            self.calc_func = self.args["calc_func"]
            self.calc_args = self.args["calc_args"]
        else:
            print("""
            ###############################################################################
            ###---   No calculator found for first image: Running setup to get it    ---###
            ###############################################################################
            """)

            self.images[0].setup()
            self.calc_func = self.images[0].calc_func
            self.calc_args = self.images[0].calc_args
            if hasattr(self.images[0], "calc_func") and hasattr(self.images[0],"calc_args"):
                print(f">>> SUCCESS: got the calculator successfully <<<\n >>> calc_func= {self.images[0].calc_func}\n >>> calc_args= {self.images[0].calc_args}")
            else:
                print("""
            ######################################################################################
            ###---   WARNING: The first image does not have a calculator defined. Leaving   ---###
            ######################################################################################
                """)
                exit(1)

    def setup(self):
        # the input and output directories find all the necessary information.        out
        self.cwd = os.getcwd()

        self.get_all_images()

        out_dir = self.args["output_directory"]

        if not os.path.exists( out_dir):
            os.mkdir(out_dir )

        self.cwd = os.getcwd()
        os.chdir(out_dir)


    def run(self):
        optimizer = BFGS(self.neb, trajectory='neb_climb.traj')
        optimizer.run(fmax=0.04)


    def get_all_images(self):

        self.get_calc()

        if self.utils.check_key(self.args, "read_traj"):
            # Read from the trajectory fild
            n_images = self.args["n_images"]
            traj_file = self.args["read_traj"]
            self.neb_images = read(f'{traj_file}@-{n_images}:')
            self.neb = NEB(self.neb_images, climb=self.climb)
        else:
            initial = self.images[ 0].structure
            final   = self.images[-1].structure

            print(f"> Initial structure:\n > {initial.positions}\n > {initial.cell}\n > {initial.symbols} ")
            print(f"> Final   structure:\n > {final.positions}\n > {final.cell}\n > {final.symbols} ")

            if len(self.images) == 2:
                print("""
                ###################################################################################
                ###---   Only 2 images: assuming initial and final:  --> Interpolating <--   ---###
                ###################################################################################
                """)


                if self.utils.check_key(self.args, "n_images"):
                    self.neb_images = [initial.copy() for _ in range(self.args["n_images"] -1)] + [final]
                else:
                    print("""
                    #####################################################################################################
                    ###---   WARNING: number of images for NEB calculation not set. Please add for interpolation   ---###
                    #####################################################################################################
                    """)
                    exit(1)
            else:
                # Get the neb images from the images object
                self.neb_images = [ image.structure for image in self.images ]


            for i,n in enumerate(self.neb_images):
                #if self.images[0].name == "GapCalc":
                # self.calc_args["directory"] = f"{i:02d}"
                # if not os.path.exists(self.calc_args["directory"]):
                #     os.mkdir(self.calc_args["directory"])
                #     self.images[0].copy_potential(self.calc_args["directory"])
                # else:
                #     self.calc_args.pop('dictionary', None)
                #elif self.calc_args.has_key("directory"):

                n.calc =  self.calc_func( **self.calc_args )

            print(f"> NEB Images ", self.neb_images)
            self.neb = NEB(self.neb_images, climb=self.climb)

            if len(self.images) == 2:
                self.neb.interpolate()
                print(">>> Interpolating...")
            print(self.neb)


    def check_convergence(self):

        return

    def get_data(self):
        pass



if __name__ == "__main__":


    test = "Vasp"
    #    test="Gap"


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


        from vasp_calculation import VaspCalc
        calculation_method = VaspCalc

        istructure = read(f"{input_directory}/POSCAR_initial", format="vasp")
        fstructure = read(f"{input_directory}/POSCAR_final", format="vasp")



        system = "CBr"

        args ={ "binary"              : binary,
                "potential_directory" : potential_directory,
                "input_directory"     : input_directory,
                "output_directory"    : output_directory,
                "input_args"          : vasp_input_args,
                "structure"           : istructure,
                "ncores"              : ncores,
                "system"              : system,
                "driver_args"         : driver_args
        }

        c1 = CalculationContainer(calculation_method,  args )

        args2 = copy.copy(args)
        args2["structure"] = fstructure
        c2 = CalculationContainer(calculation_method,  args2 )

        neb_args = {"images": [c1.method, c2.method],
                    "n_images" : 5,
                    "climb":False,
                    "input_directory"     : input_directory,
                    "output_directory"    : output_directory
                    }

        n = CalculationContainer(NEB_interface, neb_args )
        n.run()

    if test == "Gap":
        # Do neb calculation using gap

        input_directory = "input_dir"
        output_directory = "output_dir"
        potential_directory = "input_dir/gap_files"

        binary = "turbogap"
        ncores = 128

        istructure = read(f"{input_directory}/POSCAR_initial", format="vasp")
        fstructure = read(f"{input_directory}/POSCAR_final", format="vasp")

        gap_input_args = {"quip" : True}

        system = "CBr"

        args ={ "binary"              : binary,
                "potential_directory" : potential_directory,
                "input_directory"     : input_directory,
                "output_directory"    : output_directory,
                "input_args"          : gap_input_args,
                "structure"           : istructure,
                "ncores"              : ncores,
                "system"              : system
        }

        from gap_calculation import GapCalc

        calculation_method = GapCalc

        c1 = CalculationContainer(calculation_method,  args )

        args2 = copy.copy(args)
        args2["structure"] = fstructure
        c2 = CalculationContainer(calculation_method,  args2 )

        neb_args = {"images": [c1.method, c2.method],
                    "n_images" : 5,
                    "climb":False,
                    "input_directory"     : "./",
                    "output_directory"    : "neb_calc"
                    }

        n = CalculationContainer(NEB_interface, neb_args )
        n.run()
