#!/usr/bin/python3
from ase import Atoms
from ase.io import read, write
from ase.calculators.vasp import Vasp
import os, shutil, subprocess
from quippy.potential import Potential
from utils import Utils
from calculations import Calculation, CalculationUtils, CalculationContainer


class GapCalc( Calculation ):
    def __init__(self, args):
        self.name = "GapCalc"
        self.args = args # This is the calculation data
        self.utils = Utils()
        self.calc_utils = CalculationUtils()

        self.result = {}
        self.structure = self.args.structure

    def __str__(self):
        return f"{self.name} = (args = {self.args.__str__()}, result = {self.result})"



    def find_potential_file(self, directory):
        # First check this current directory
        pot_file = self.utils.check_file_dir_subdir(f'{self.args.system}.xml', dir='.', subdir="gap_files")

        if pot_file is None:
            pot_file = self.utils.check_file_dir_subdir(f'{self.args.system}.xml', dir=directory, subdir="gap_files")

        if pot_file is None:
            print(f"""
            ###################################################################################################################################################
            ###---   FATAL: Could not find the potential file {self.args['system']}.xml in ./ or ./gap_files or {directory} or {directory}/gap_files   ---###
            ###################################################################################################################################################

            """)
            exit(1)
        else:
            return os.path.abspath(pot_file)

    def copy_potential(self, dir):
        out_dir = self.args.output_directory
        pot_file = self.find_potential_file(out_dir)
        name = pot_file.split("/")[-1]
        shutil.copy(pot_file, f"{dir}/{name}")

    def setup(self):
        # Copy gap files from directory to where the calculation is

        out_dir = self.args.output_directory
        pot_file = self.find_potential_file(out_dir)

        self.path = os.path.abspath(out_dir)
        self.pot_path = os.path.abspath(pot_file)

        if self.utils.check_key(self.args.input_args, "quip"):
            self.setup_quip(pot_file, out_dir)
        else:
            self.setup_turbogap()



    def run(self):
        # Create the input file for the directory and then compute
        self.result = self.calc_utils.run(self.structure, self.args)

        return self.result


    def setup_quip(self, pot_file, out_dir):
        print(">>>   Setting up QUIP gap calculation <<<")
        gap = Potential(param_filename=pot_file, directory = out_dir)
        self.structure.set_calculator(gap)
        print(">>>   Assigning calculator <<<")
        self.calc_func = Potential
        self.calc_args = {"param_filename":self.pot_path, "directory" : out_dir}



    def setup_turbogap(self):
        pass

    def get_data(self):
        pass

    def save(self, prefix=''):
        # Write the json file
        from ase.io import jsonio
        dct = self.structure.calc.results  # Get the calculator in a dictionary format
        dct_extra = self.structure.calc.extra_results  # Get the calculator in a dictionary format

        if len(prefix) == 0:
            prefix = self.name
        name = self.utils.get_save_name(self.path, self.result, prefix)
        jsonio.write_json(name, dct)
        jsonio.write_json(name.replace( ".json", "_extra.json" ), dct_extra)


    def save_get_forces(self, atoms, name="", dir="00"):


        # The below is modified from ASE, this is so we can get all the OUTCAR files saved somewhere during calculation
        def get_forces(apply_constraint=True, md=False):
            """Calculate atomic forces.

            Ask the attached calculator to calculate the forces and apply
            constraints.  Use *apply_constraint=False* to get the raw
            forces.

            For molecular dynamics (md=True) we don't apply the constraint
            to the forces but to the momenta. When holonomic constraints for
            rigid linear triatomic molecules are present, ask the constraints
            to redistribute the forces within each triple defined in the
            constraints (required for molecular dynamics with this type of
            constraints)."""

            if atoms._calc is None:
                raise RuntimeError('Atoms object has no calculator.')
            forces = atoms._calc.get_forces(atoms)

            if apply_constraint:
                # We need a special md flag here because for MD we want
                # to skip real constraints but include special "constraints"
                # Like Hookean.
                for constraint in atoms.constraints:
                    if md and hasattr(constraint, 'redistribute_forces_md'):
                        constraint.redistribute_forces_md(atoms, forces)
                    if not md or hasattr(constraint, 'adjust_potential_energy'):
                        constraint.adjust_forces(atoms, forces)

            #self.structure.calc.write_json(name)
            from ase.io import jsonio
            dct = atoms.calc.results  # Get the calculator in a dictionary format
            dct_extra = atoms.calc.extra_results  # Get the calculator in a dictionary format

            prefix=f"{self.name}_calc"
            if len(prefix) == 0:
                prefix = self.name
            name = self.utils.get_save_name(self.path, self.result, prefix)
            jsonio.write_json(name, dct)
            jsonio.write_json(name.replace( ".json", "_extra.json" ), dct_extra)


            return forces


        return get_forces
