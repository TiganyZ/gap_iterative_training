#!/usr/bin/python3
from ase import Atoms
from ase.io import read, write
from ase.calculators.vasp import Vasp
import os, shutil, subprocess
from utils import Utils
from calculations import Calculation, CalculationUtils, CalculationContainer


class VaspCalc( Calculation ):

    def __init__(self, args):
        self.name = "VaspCalc"
        self.args = args
        self.utils = Utils()
        self.calc_utils = CalculationUtils()
        self.result = {}
        self.counter = 0

        self.structure = self.args.structure

    def check_vasp_input_args(self):

        default_input_args = {"prec" : "Accurate",
                              "ibrion" : 2,
                              "encut" : 650,
                              "ediff" : 1.0e-08,
                              "ediffg" : -1e-03,
                              "ismear" : 0,
                              "sigma" : 0.1,
                              "algo" : "fast",
                              "lwave" : False,
                              "lcharg" : False,
                              "istart" : 0,
                              "nsw" : 0,
                              "ncore" : 16 ,
                              "isif" : 2,
                              "kspacing" : 0.25,
                              "kgamma" : True,
                              }

        for k,v in self.args.input_args.items():
            if k in default_input_args.keys():
                # Compare the values
                default = default_input_args[k]

                if default /= v:
                    print(f"""
                    ########################################################################################################################
                    ###---   VASP KEY CHECK: Input for {k} = {v} is not equal to default {k} = {default}: Make sure you know why!!!   ---###
                    ########################################################################################################################
                    """)
            else:
               print(f"""
               ####################################################################################################################
               ###---   VASP KEY CHECK: {k} is not in default arg list for standard calculation. Make sure this is correct   ---###
               ####################################################################################################################
               """)

        
    def __str__(self):
        return f"{self.name} = (args = {self.args.__str__()}, result = {self.result})"

    def setup(self):
        # Copy gap files from directory to where the calculation is
        out_dir = self.args.output_directory
        pot_dir = self.args.potential_directory
        self.pot_path = os.path.abspath(f"{pot_dir}/POTCAR" )


        if self.args.make_dirs:
            self.copy_potential(out_dir)
            self.path = os.path.abspath(out_dir)
        else:
            self.path = os.path.abspath("./")


        self.check_vasp_input_args()

        self.args.input_args["directory"] = self.path
        self.structure.calc = Vasp( **self.args.input_args )
        self.calc_func = Vasp
        self.calc_args = self.args.input_args

        self.create_run_environment(out_dir)


    def copy_potential(self, dir):
        out_dir = self.args.output_directory
        pot_dir = self.args.potential_directory
        self.utils.check_copy_file(pot_dir, "POTCAR", dir)

    def create_run_environment(self, out_dir):
        self.args.driver_args["ncores"] = self.args.binary
        self.args.driver_args["binary"] = self.args.ncores

        if self.args.batch:
            driver_template = self.utils.get_driver_template(self.args.driver_args)

        else:
            binary = self.args.binary
            ncores = self.args.ncores
            # Make the run_vasp for the number of cores that we want
            with open(f"{out_dir}/run_vasp.py", 'w') as f:
                f.write(f"""
import os
exitcode = os.system('srun -n {ncores} {binary}')
    """)

            cwd = os.getcwd()

            os.environ["VASP_SCRIPT"]=os.path.abspath(f"{out_dir}/run_vasp.py")





    def run(self):


        # Create the input file for the directory and then compute
        # self.utils.check_keys(self.args, keys=( "structure", "input_args" ) )

        # self.structure.calc = Vasp( **self.args.input_args )
        # self.calc = Vasp( **self.args.input_args )

        # print(self.structure.calc)
        self.result = self.calc_utils.run(self.structure, self.args)

        # print("VASP ENERGY RESULT: ", self.result)
        # self.result = {"result":self.result}#self.calc_utils.run(self.structure, self.args)

        return self.result



    def save(self, prefix=""):

        if len(prefix) == 0:
            prefix = self.name
        name = self.utils.get_save_name(self.path, self.result, prefix)
        self.structure.calc.write_json(name)


    def check_convergence(self):
        if os.path.exists(f"{self.path}/OUTCAR"):
            self.converged = self.structure.calc.read_convergence()
            print(f">>> CHECK CONVERGENCE {self.name}: Has calculation converged? {self.converged}")

        else:
            print("""
            ###################################################################################################
            ###---   WARNING: Convergence for the calculation has not been reached. Check OUTCAR file.   ---###
            ###################################################################################################

""")


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

            filename = self.utils.get_save_name(f"{dir}/jsons", {}, f"{self.name}_calc")
            atoms.calc.write_json(f"jsons/{filename}")
            #self.utils.save_file_in_dir(filename, dir, "jsons" )

            filename = self.utils.get_save_name(f"{dir}/outcars", {}, "OUTCAR", ext="")
            shutil.copy(f"{dir}/OUTCAR", f"{dir}/outcars/{filename}" )
            #self.utils.save_file_in_dir("OUTCAR", dir, "outcars" )

            filename = self.utils.get_save_name(f"{dir}/images", {}, f"{self.name}_{dir}", ext=".xyz")
            write( f"{dir}/images/{filename}", atoms, format="extxyz" )
                #self.utils.save_file_in_dir(filename, dir, "images" )


            return forces

        self.counter += 1
        return get_forces

    def get_data(self):
        pass
