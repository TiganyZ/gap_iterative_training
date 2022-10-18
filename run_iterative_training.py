#!/usr/bin/python3
"""This is the main file for the running of iterative training.

The purpose of this program of to automate the iterative training of a
GAP using ase. This general purpose use will then be used to
iteratively train a NEB.

Can flesh out some bits here


"""
#!/usr/bin/python3

from abc import ABC, abstractmethod
from typing import Type, TypeVar
from datetime import datetime
from ase import Atoms
from ase.io import read, write
from ase.calculators.vasp import Vasp
import os, shutil, subprocess

from utils import Utils
from process_calculation import GAP_to_VASP, VASP_to_GAP











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
