#!/usr/bin/python3

from calculations import CalculationData, CalculationContainer
from neb_interface import NEB_data, NEB_interface
import os, shutil, subprocess, pickle, sys
from ase.io import read, write


# These are actually the results of a relaxed gap calculation
if len(sys.argv) >= 3:
    files = [sys.argv[1], sys.argv[2]]
else:
    print("!!! WARNING: There have two files specified for the NEB calculation !!!")
    exit(1)
    


head, tail = os.path.split(files[0])
base = tail.split('.')[0]

istructure = read(files[0])
fstructure = read(files[1])


gap_input_args = {"quip" : True,
                  "pot1" : "carbon.xml",
                  "pot2" : "CBr.xml"}


args = CalculationData( binary              = "-",
                        input_directory     = "./",
                        potential_directory = "./gap_files",
                        output_directory    = f"output_NEB_{base}",
                        structure           = istructure,
                        input_args          = gap_input_args,
                        ncores              = 1,
                        system              = "C+Br", 
                        run_calc_type       = "energy"
                       )

# Do neb calculation using gap


from gap_calculation import GapCalc

calculation_method = GapCalc

c1 = CalculationContainer(calculation_method,  args )



args = CalculationData( binary              = "-",
                        input_directory     = "./",
                        potential_directory = "./gap_files",
                        output_directory    = f"output_NEB_{base}",
                        structure           = fstructure,
                        input_args          = gap_input_args,
                        ncores              = 1,
                        system              = "C+Br", 
                        run_calc_type       = "energy"
                       )

c2 = CalculationContainer(calculation_method,  args )

neb_args = NEB_data(images = [c1.method, c2.method],
                    n_images = 5,
                    climb = False
                    )

with open("neb_args.pkl", "wb") as f:
    pickle.dump( neb_args, f)

n = CalculationContainer(NEB_interface, neb_args )
n.method.name = f"{n.method.name}_{base}"
n.run()

# Running the NEB calculation there will be a neb_calc folder by default

relfiles = []
for i in range(1, neb_args.n_images - 1):
    # copy the images to the current folder
    idir=f"0{i}/images"
    # Now get the latest file

    if not os.path.exists(idir):
        os.makedirs(idir, exist_ok=True)

    l = (os.listdir(idir))
    counter =0
    lm = l[0]
    nlp=0
    for j, lj in enumerate(l):
        try:
            nl = int( (lj.split('.')[0]).split('_')[-1] )
            if nl > nlp:
                lm = lj
                nlp = nl
            else:
                continue
        except ValueError:
            continue
    print(f">> Last image found: {lm} in {idir}\n>> Copying to ./relaxed_neb_image_0{i}_{base}.xyz")
    
    final_relaxed = lm #max(, key=os.path.getctime)
    
    shutil.copy(f"{idir}/{final_relaxed}", f"../relaxed_neb_image_0{i}_{base}.xyz")
    # Add this to the list of vasp calculations to be done.
    # Now pass on to vasp, if there are two potentials, then make two vasp calculations
    # Delete the br species
    relfiles.append(f"relaxed_neb_image_0{i}_{base}.xyz")
    
    element = 'Br'
    atoms = read(f"../relaxed_neb_image_0{i}_{base}.xyz")
    del atoms[[atom.index for atom in atoms if atom.symbol == element]]
    write(f"../relaxed_neb_image_0{i}_{base}_no_{element}.xyz", atoms)
    relfiles.append(f"relaxed_neb_image_0{i}_{base}_no_{element}.xyz")


os.chdir("../")
vfile = "run_vasp_of_gap.py"
for rstr in relfiles:
    output = "vasp_calculations_to_run.sh"
    command = f"python3.9 {vfile} {rstr}\n"

    if os.path.exists(output):
        mode = 'a'
    else:
        mode = 'w'

    with open(output, mode) as f:
        f.write( command )

    
    
