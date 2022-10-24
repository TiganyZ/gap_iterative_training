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
    def __init__(self, neb_args):

        self.name = "NebCalc"
        self.images = self.neb_args["images"]
        self.neb_args = neb_args
        self.utils = Utils()

        if self.utils.check_key(neb_args, "climb"):
            self.climb = self.neb_args["climb"]
        else:
            self.climb = False


    def get_calc(self):
        if hasattr(self.images[0].structure, "calc"):
            self.calc = self.images[0].structure.calc
        else:
            print("""
            ###############################################################################
            ###---   No calculator found for first image: Running setup to get it    ---###
            ###############################################################################
            """)

            self.images[0].setup()
            os.chdir("../")
            if hasattr(images[0].structure, "calc"):
                print(">>> SUCCESS: got the calculator successfully <<<")
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

        if self.utils.check_key(self.neb_args, "read_traj"):
            # Read from the trajectory fild
            n_images = self.neb_args["n_images"]
            traj_file = self.neb_args["read_traj"]
            self.neb_images = read(f'{traj_file}@-{n_images}:')
            self.neb = NEB(self.neb_images, climb=self.climb)
        else:
            initial = self.images[ 0].structure
            final   = self.images[-1].structure

            if len(self.images) == 2:
                print("""
                ###################################################################################
                ###---   Only 2 images: assuming initial and final:  --> Interpolating <--   ---###
                ###################################################################################
                """)


                if self.utils.check_key(self.neb_args, "n_images"):
                    self.neb_images = [initial.copy() for i in range(self.neb_args["n_images"] -1)] + [final]
                    self.neb_images = [n.set_calculator(self.calc) for n in self.neb_images]
                    self.neb = NEB(self.neb_images, climb=self.climb)
                    self.neb.interpolate()
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
                self.neb_images = [n.set_calculator(self.calc) for n in self.neb_images]
                self.neb = NEB(self.neb_images, climb=self.climb)


    def get_data(self):
        pass



if __name__ == "__main__":


    test = "Vasp"
    test="Gap"



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
                    "input_directory"     : input_directory,
                    "output_directory"    : output_directory
                    }

        n = CalculationContainer(NEB_interface, neb_args )
        n.run()