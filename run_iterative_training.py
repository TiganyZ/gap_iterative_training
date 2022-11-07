#!/usr/bin/python3
from calculations import CalculationData, CalculationContainer
from train_gap import Train
import os, shutil, subprocess
from ase.io import read, write


name="graphene_relax"

system = "CBr"

d="outcars_combined"
outcars_combined = [ os.path.join(d,i) for i in os.listdir(d) ]
d="outcars_isolated"
outcars_isolated = [ os.path.join(d,i) for i in os.listdir(d) ]

outcars= {"combined": outcars_combined,
          "isolated": outcars_isolated }

previous_gap = "~/Documents/Br_C_GAP/for_tigany/dimer/refit_dimer_to_spin-orbit/CBr/so_z_antiparallel_new_8/gap/gap_files"  #"../../../prev_gap/gap_files"

info = {"previous_database": "./train_NPC_graphite_dimer_so_z.xyz",

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
        "numbers" : { "C" : 6,             "Br" : 35       },
        "masses"    : { "C" : 12.01,         "Br" : 79.90412 },
        "e0"      : { "C" : -.16138053,    "Br" : 0.0      },
        "l_max" : 8,
        "n_max" : [8, 8]

        }


t = Train( name,
           system,
           outcars,
           info,
           delta_gap = True,
           previous_gap = previous_gap
)

t.setup()

t.run()
