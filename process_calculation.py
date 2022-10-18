from abc import ABC, abstractmethod
from typing import Type, TypeVar
from datetime import datetime
from ase import Atoms
from ase.io import read, write
from ase.calculators.vasp import Vasp
import os, shutil, subprocess

from calculations import Calculation, GapCalc, VaspCalc, CalculationContainer


class ProcessCalculation(ABC):

    @abstractmethod
    def extract_data():
        pass

    @abstractmethod
    def create_converted_calculation():
        pass



class GAP_to_VASP(ProcessCalculation):

    def __init__(self, calc: Type[CalculationContainer], args: dict):
      self.calc = calc
      self.args = args

      self.structure = self.calc.args["structure"]

      if self.calc.run_energy:
          # Then extract the data






class VASP_to_GAP(ProcessCalculation):

    def __init__(self, calc: Type[CalculationContainer], args: dict):
      self.calc = calc
      self.args = args
