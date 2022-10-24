#!/usr/bin/python3

from abc import ABC, abstractmethod
from typing import Type, TypeVar
from datetime import datetime
from ase import Atoms
from ase.io import read, write
from ase.calculators.vasp import Vasp
import os, shutil, subprocess, re
import numpy as np
from utils import Utils


class Train:
    # Not creating abstract class because we will only be training GAPs

    def __init__(self, system, outcars, info, sigma_e = 0.0002, sigma_v = 0.02, delta_gap = False, name="", train_file_name="train.xyz", verbosity=100):
        self.system = system
        self.outcars = outcars # Dictionary of outcars { "C+element": [out1, out2], "C":[out1, out2] }
        self.delta_gap = delta_gap
        self.name = name
        self.elements = re.findall('[A-Z][^A-Z]*', system)
        self.train_file_name = train_file_name
        self.verbosity = verbosity

        self.info = info # This has information on things like the path of the current training base, the parameters (sigma_e and sigma_v) and names of the databases

        self.dbs_name = f"{'_'.join(self.elements)}" + self.name

        self.masses = ' '.join([ str((self.info["mass"])[e]) for e in self.elements ])
        self.e0     = ' '.join([ str((self.info["e0"])[e]) for e in self.elements ])
        self.number = ' '.join([ str((self.info["number"])[e]) for e in self.elements ])


        self.utils = Utils()


    def create_dbs(self):
        append = False
        keys = tuple(self.outcars.keys())
        for out_files in zip(*tuple(self.outcars.values())):

            if np.all([ os.path.isfile(file) for file in out_files ]):

                if self.delta_gap:
                    atoms1 = read( out_files[ keys.index("combined") ] )
                    atoms2 = read( out_files[ keys.index("isolated") ] )

                    atoms = self.get_dbs_data_delta_gap(atoms1, atoms2)
                else:
                    atoms1 = read( out_files["combined"] )
                    atoms = self.get_dbs_data_gap(atoms1)

                write(f"{self.dbs_name}.xyz" , atoms, append=append)
                append = True

        # Now append this to a previous file
        self.utils.check_key(self.info, "previous_database", stop=True)

        with open(self.train_file_name, "w") as f:
            out = subprocess.run(["cat", f"{self.info['previous_database']}", f"{self.dbs_name}.xyz"], stdout=f)
            self.utils.check_subprocess(out)

        self.get_data_from_database(self.train_file_name)



    def get_data_from_database(self, dbs_file):

        commands = f"grep config_type {dbs_file} | awk -F\"config_type\" '/config_type/{{print $2}}' | sed 's/=//g' | awk '{{print$2}}' "
        out = self.utils.piped_subprocess(commands)

        # out has all the filenames of the config types
        self.config_types = out.split()

        self.config_dict = {}
        checked = []
        for config in self.config_types:
            if config in checked:
                self.config_dict[config] += 1
            else:
                self.config_dict[config]  = 1
                checked.append(config)

        if self.verbosity > 50:
            print(f" >>> Get data from database {dbs_file} <<<\n  > Synopsis\n  {config_dict}")

        return self.config_dict

    def get_dbs_data_gap(self, atoms1):

        pos = atoms1.get_positions()
        cell = atoms1.get_cell()
        symb = atoms1.get_chemical_symbols()
        n1 = len(atoms1)
        forces = np.zeros([n1, 3])
        forces[0:n1] = atoms1.get_forces()[0:n1]
        e = atoms1.get_potential_energy(force_consistent=True)
        v = atoms2.get_volume()
        virial = -v*(atoms1.get_stress(voigt=False))
        atoms = Atoms(symb, cell=cell, positions=pos, pbc=True)
        atoms.info["free_energy"] = e
        atoms.info["virial"] = virial
        atoms.info["config_type"] = self.dbs_name
        atoms.set_array("forces", forces)
        return atoms


    def get_dbs_data_delta_gap(self, atoms1, atoms2):
        pos = atoms2.get_positions()
        cell = atoms2.get_cell()
        symb = atoms2.get_chemical_symbols()
        n1 = len(atoms1)
        n2 = len(atoms2)
        forces = np.zeros([n1+n2, 3])
        forces[0:n1] = atoms2.get_forces()[0:n1] - atoms1.get_forces()
        forces[n1:n2] = atoms2.get_forces()[n1:n2]
        e = atoms2.get_potential_energy(force_consistent=True) - atoms1.get_potential_energy(force_consistent=True)
        v = atoms2.get_volume()
        virial = -v*(atoms2.get_stress(voigt=False) - atoms1.get_stress(voigt=False))
        atoms = Atoms(symb, cell=cell, positions=pos, pbc=True)
        atoms.info["free_energy"] = e
        atoms.info["virial"] = virial
        atoms.info["virial"] = virial
        atoms.info["config_type"] = self.dbs_name
        atoms.set_array("forces", forces)
        return atoms



    def add_tags(self):
        # This adds the structures to the add tags template
        # Need the names of the other structures and then the additional errors
        se = self.info["sigma_e"]
        sv = self.info["sigma_v"]

        new_struc = f"{self.dbs_name}"
        if (new_struc in se.keys()) or (new_struc in sv.keys()):
            print("""
            ##########################################################################################################
            ###---   WARNING!! The name of the additional structures is already in the database!! Change this   ---###
            ##########################################################################################################
            """)
        else:
            se[f"{self.dbs_name}"] = sigma_e
            sv[f"{self.dbs_name}"] = sigma_v


        input_dict = { "sigma_e" : se, "sigma_v" : sv }
        with open("add_tags_template.py", 'r') as f:
            file_with_data = f.read() % input_dict

        with open("add_tags_new.py", "w") as f:
            f.write(file_with_data)

        out = subprocess.run( ["python3", "add_tags_new.py"] )
        self.utils.check_subprocess(out)


    def make_gap_files(self):
        out = subprocess.run( ["python3", "make_gap_files.py", f"{self.system}.xml", f"{self.system}.gap"] )
        self.utils.check_subprocess(out)



    def create_compress_data(self):

        with open("compress_indices_temp.py", 'w') as f:
            commands = f"cat compress_indices.py | sed 's/nmax = placeholder/nmax = {self.info['n_max']} /g | sed 's/lmax = placeholder/lmax = {self.info['l_max']} /g'  "
            out = self.utils.piped_subprocess(commands, file = f )

        cmd = "python3 compress_indices_temp.py"
        out = subprocess.run( cmd.split )
        self.utils.check_subprocess(out)



    def create_train_sh(self):

        e0 = { el:self.e0[i] for i,el in enumerate(self.element) }

        gap_2b_terms = ""
        checked = []
        for i,el in enumerate(self.elements):
            for j,el2 in enumerate(self.elements):
                if i == 0 and j == 0:
                    gap_2b_terms += "distance_2b Z1={self.numbers[i]} Z2={self.numbers[j]} cutoff=5.5 n_sparse=40 covariance_type=ard_se delta=0.5 theta_uniform=0.5 sparse_method=uniform add_species=F "
                elif i == len(elements)-1 and j == len(elements)-1:
                    break
                else:
                    if el2+el not in checked:
                        gap_2b_terms += ": \ \n distance_2b Z1={self.numbers[i]} Z2={self.numbers[j]} cutoff=5.5 n_sparse=40 covariance_type=ard_se delta=0.5 theta_uniform=0.5 sparse_method=uniform add_species=F "
                        checked.append(el+el2)

        if len(elements) > 1:
            central_index = 2
        else:
            central_index = 1


        with open("train.sh", "w") as f:
            f.write(f"""
#!/bin/bash

energy="free_energy"
forces="forces"
virial="virial"


gap_fit atoms_filename=train_tagged.xyz \
        gap = {{ {gap_2b_terms}
            soap_turbo l_max={self.info['l_max']} alpha_max={{{self.info['n_max'][0]} {self.info['n_max'][1]}}} atom_sigma_r={{{0.5 0.5}}} atom_sigma_t={{{0.5 0.5}}} \
                    atom_sigma_r_scaling={{{0. 0.}}} atom_sigma_t_scaling={{{0. 0.}}} zeta=4 rcut_soft=5.0 rcut_hard=5.5 \
                    basis=poly3gauss scaling_mode=polynomial amplitude_scaling={{{1.0 1.0}}} n_species=2 \
                    species_Z={{{self.numbers[1]} {self.numbers[0]}}} radial_enhancement={{{1}}} compress_file=compress.dat central_index={central_index} \
                    central_weight={{{1.0 1.0}}} add_species=F \
                    config_type_n_sparse={self.config_dict} \
                    delta=0.1 f0=0.0 covariance_type=dot_product sparse_method=cur_points }} \
                 default_sigma={{0.001 0.1 0.1 0.1}} energy_parameter_name=$energy force_parameter_name=$forces \
                 force_mask_parameter_name=force_mask virial_parameter_name=$virial sparse_jitter=1.0e-8 do_copy_at_file=F \
                 sparse_separate_file=T gp_file={self.system}.xml openmp_chunk_size=100 e0={e0}
            """)

        return

    def create_turbogap_input():

        with open("input", 'w') as f:
            f.write(f"""
! Species-specific info
atoms_file = 'train.xyz'
pot_file = 'gap_files/{self.system}.gap'
n_species = {len(self.elements)}
species = {' '.join(self.elements)}
masses = {self.masses}
e0 = {self.e0}
""")


    def wrap_function(self, function, message):
        print(f"> Train.setup:          {message} ")
        function()
        print(f"> Train.setup: SUCCESS: {message} \n")

        return




    def setup(self):

        self.wrap_function( self.create_dbs,      "Creating database")
        self.wrap_function( self.add_tags,        "Adding tags")
        self.wrap_function( self.make_gap_files,  "Converting gap_fit to turbogap")
        self.wrap_function( self.create_train_sh, "Creating train.sh file")




    def create_workflow(self, args):
        # This is what is usually the workflow shell script
        return


if __name__ == "__main__":


    system = "CBr"

    outcars= {"combined": ["OUTCAR_1", "OUTCAR_2"]}

    info = {"previous_database": "./train.xyz",
            "sigma_e" : {"default": 0.001,
                         "nanoporous": 0.002,
                         "graphite_v1": 0.0005,
                         "graphite_v2": 0.0005,
                         "graphite_v3": 0.0005,
                         "graphite_v4": 0.0005,
                         "graphite_v5": 0.0005,
                         "graphite_v6": 0.0005,
                         "graphite_v7": 0.0005,
                         "dimer": 0.0005},
            "sigma_v" : {"default": 0.1,
                       "nanoporous": 0.2,
                       "graphite_v1": 0.05,
                       "graphite_v2": 0.05,
                       "graphite_v3": 0.05,
                       "graphite_v4": 0.05,
                       "graphite_v5": 0.05,
                       "graphite_v6": 0.05,
                       "graphite_v7": 0.05,
                         "dimer": 0.05},
            "numbers" : { "C" : 12,  "Br" : 35 },
            "masses"  : { "C" : 12,  "Br" : 35 },
            "e0"      : { "C" : 0.0, "Br" : 0.0 },
            "l_max" : 8,
            "n_max" : [8, 8]

            }


    t = Train( system, outcars, info )
