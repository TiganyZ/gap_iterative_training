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
        self.name = f"{args.system}_GapCalc"
        self.args = args # This is the calculation data
        self.utils = Utils()
        self.calc_utils = CalculationUtils()

        self.result = {}
        self.structure = self.args.structure

    def __str__(self):
        return f"{self.name} = (args = {self.args.__str__()}, result = {self.result})"



    def find_potential_file(self, file, directory):
        # First check this current directory
        pot_file = self.utils.check_file_dir_subdir(file, dir='.', subdir="gap_files")

        if pot_file is None:
            pot_file = self.utils.check_file_dir_subdir(file, dir=directory, subdir="gap_files")

        if pot_file is None:
            print(f"""
            ###################################################################################################################################################
            ###---   FATAL: Could not find the potential file {file} in ./ or ./gap_files or {directory} or {directory}/gap_files   ---###
            ###################################################################################################################################################

            """)
            exit(1)
        else:
            return os.path.abspath(pot_file)


    def get_potential_files(self, out_dir):
        self.pot_files = []
        for key in ["pot1", "pot2"]:
            if self.utils.check_key(self.args.input_args, key):
                self.pot_files.append( self.find_potential_file(self.args.input_args[key], out_dir))
                self.args.input_args[key] = self.pot_files[-1]
            else:
                self.pot_files.append( self.find_potential_file(f"{self.system}.xml", out_dir))



    def copy_potential(self, dir):
        out_dir = self.args.output_directory

        for pot_file in self.pot_files:
            name = pot_file.split("/")[-1]
            shutil.copy(pot_file, f"{dir}/{name}")


    def setup(self):
        # Copy gap files from directory to where the calculation is

        out_dir = self.args.output_directory
        self.pot_path = os.path.abspath(self.args.potential_directory)
        self.path = os.path.abspath(out_dir)
        #        self.pot_path = os.path.abspath(pot_file)

        if self.utils.check_key(self.args.input_args, "quip"):
            self.args.input_args.pop("quip")
            self.setup_quip(out_dir)

        else:
            self.setup_turbogap()



    def run(self):
        # Create the input file for the directory and then compute
        self.result = self.calc_utils.run(self.structure, self.args, name=self.name, path=self.path)

        return self.result


    def setup_quip(self, out_dir):
        print(">>>   Setting up QUIP gap calculation <<<")
        self.get_potential_files(out_dir)

        if len(self.pot_files) == 2:
            print(self.pot_files)
            pot1 = Potential( param_filename = self.pot_files[0])
            pot2 = Potential( param_filename = self.pot_files[1])

            gap = Potential(args_str = "Sum", pot1=pot1, pot2=pot2 )
            self.calc_args = { "args_str" : "Sum",
                               "pot1" : pot1,
                               "pot2" : pot2,
                               "directory" : self.path}
        else:
            gap = Potential( param_filename=self.pot_file[0]) #, directory = out_dir)
            self.calc_args = {"param_filename":self.pot_file[0], "directory" : self.path}
        self.structure.set_calculator(gap)
        print(">>>   Assigning calculator <<<")
        self.calc_func = Potential




    def setup_turbogap(self):
        pass

    def get_data(self):
        pass

    def save_gap_files(self, atoms, dir="."):
        from ase.io import jsonio
        dct = self.structure.calc.results  # Get the calculator in a dictionary format
        dct_extra = self.structure.calc.extra_results  # Get the calculator in a dictionary format

        prefix = self.name


        name =self.utils.get_save_name(f"{dir}/jsons", self.result, prefix)
        jsonio.write_json(os.path.join(f"{dir}/jsons", name), dct)
        jsonio.write_json(os.path.join(f"{dir}/jsons", name.replace( ".json", "_extra.json" )), dct_extra)

        filename = self.utils.get_save_name(f"{dir}/images", {}, f"{prefix}", ext=".xyz")
        write( f"{dir}/images/{filename}", atoms, format="extxyz" )

    
    def save(self, prefix=''):
        # Write the json file

        if hasattr(self, "result"):
            if self.utils.check_key( self.result, "optimized_structure" ):
                self.save_gap_files(self.result["optimized_structure"], dir=self.path)
                return None

        self.save_gap_files(self.structure, dir=self.path)
        return None



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


            self.save_gap_files(atoms, dir)

            # filename = self.utils.get_save_name(dir, {}, f"{self.name}_{dir}", ext=".json")
            # jsonio.write_json(f"{dir}/{filename}", dct)
            # self.utils.save_file_in_dir(filename, dir, "jsons" )

            # jsonio.write_json(f"{dir}/{filename.replace( '.json', '_extra.json' )}", dct_extra)
            # self.utils.save_file_in_dir(filename.replace( '.json', '_extra.json' ), dir, "extra_jsons" )


            # filename = self.utils.get_save_name(dir, {}, f"{self.name}_{dir}", ext=".xyz")
            # write( f"{dir}/{filename}", atoms, format="extxyz" )

            # self.utils.save_file_in_dir(filename, dir, "images" )

            return forces


        return get_forces
