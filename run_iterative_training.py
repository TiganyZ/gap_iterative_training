#!/usr/bin/python3
"""This is the main file for the running of iterative training.

The purpose of this program of to automate the iterative training of a
GAP using ase. This general purpose use will then be used to
iteratively train a NEB.

Can flesh out some bits here


"""
from abc import ABC, abstractmethod
from typing import Type, TypeVar
from datetime import datetime
from ase import Atoms
from ase.io import read, write
import os, shutil
#################################
###---   Structure Class   ---###
#################################
# > Perhaps with this, one can make it dependent on the ase
# > classes/subclasses used for the objects.


class Calculation(ABC):

    @abstractmethod
    def setup():
        pass:

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


    def setup(self):
        # Copy gap files from directory to where the calculation is
        potdir = args["potential_directory"]
        outdir = args["output_directory"]

        if os.path.exists(f"{potdir}/POTCAR"):
            shutil.copy( f"{potdir}/POTCAR", f"{outdir}/" )
        else:
            print(f"""
            ############################################################
            ###---   WARNING! There is no path {potdir/POTCAR}!   ---###
            ############################################################\n\n  Rectify this, I'm leaving you...""")
            exit(1)



        shutil.copytree( args[ "potential_directory" ], args["output_directory"] )



class GAPCalc( Calculation ):

    def __init__(self, args):
        self.name = "GAPCalc"
        self.args = args

    def setup(self):
        # Copy gap files from directory to where the calculation is
        shutil.copytree( args[ "potential_directory" ], args["output_directory"] )





class CalculationContainer:
    def __init__(self, calculation_method: Type[Calculation],
                 binary: str, potential_directory: str, output_directory: str):


        self.binary = binary
        self.potential_directory = potential_directory
        # self.images = sorted(os.listdir(self.potential_directory))

        self.output_directory = output_directory

        args = { "binary": binary, "potential_directory": potential_directory, "output_directory":output_directory }

        self.method = calculation_method(args)
        self.method_name = self.method.name

        now = datetime.now()
        self.dt = now.strftime("%Y-%m-%d--%H-%M-%S")



    def create_output_directory(self):
        dir_name = f"{self.method_name}_{self.dt}_{self.image_directory}"
        if not os.path.exists(dir_name):
            os.mkdir(dir_name)
        return dir_name


    def run(self):
        self.create_output_directory()

        self.method.setup()
        self.method.get_energy()
        self.method.get_data()


    
struc = Structure( species, positions, periodicity )
gap   = calculation( calculation_type="GAP", source="some/directory/gap_files"  )

# Now compose the structures and the calculation type
# > Run GAP calculation with structure
r = Run( gap, structure )
r.run() # This will create a directory for the files etc


# Generalise the the RunNEB class

n = RunNEB( [ Run(gap, struc1), Run(gap, struc2)... ]  )

# Now one can get the data from the NEB runs

n.get_final_images() # This has the energy and the configuration
tst = n.get_tst() # This is a structure type


# Now feed this into a DFT calculation

dft = calculation(calculation_type = "VASP")
d = Run( dft, tst )
d.run()

# Now get the data from the DFT run and the gap info to retrain
data = ProcessData( type="FitGAP", d  ) # Know that there is a DFT calculation so will extract energy from there
data.extract()

retrain = RefitGAP( data, gap )
