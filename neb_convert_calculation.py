#!/usr/bin/python3

from ase.io import read, write
import pickle, os, shutil


"""

Directory structure for NebCalc is as follows:

> results_{timestamp}
-> neb_calc
--> 01 02 03 ... 0N
    [if VaspCalc]
---> images outcars
    [if GapCalc]
---> images jsons extra_jsons

Now need to get the data for a particular run of the NEB: this can be run within the results_{timestamp} directory
This will create a directory with image files from the NEB

> results_{timestamp}

"""

# Just need to extract the structure files from the xyz in the images and then make a new object

class ConvertNEB:
    def __init__(self, results_directory, output_directory, create,  neb_dir="neb_calc"):
        self.results_directory = results_directory
        self.output_directory = output_directory
        self.neb_directory = neb_dir
        self.create = create


    def get_last_image(self, directory):
        l = sorted(os.listdir(directory))
        image = l[-1]
        print(f" Last file in directory {directory} = {image}")
        return os.join.path(directory, image)

    def get_image_dirs(self, path):
        image_dirs = []
        for l in sorted(os.listdir(path)):
            if os.isdir(l):
                try:
                    n = int(l)
                    print(f"{l} is image dir")
                    image_dirs.append(l)
                except ValueError:
                    print(f"{l} is not image dir")
        return image_dirs

    def copy_images(self, image_dirs):
        image_paths = []
        for i,id in enumerate(image_dirs):
            if not os.path.exists(self.output_directory):
                os.mkdir(self.output_directory)
            last_image = get_last_image(f"{id}/images")
            filename = f"image_{i}.xyz"
            shutil.copy(last_image, f"{self.output_directory}/{filename}")
            image_paths.append(f"{self.output_directory}/{filename}")


    def calculation_method(self):
        if self.create == "Vasp":
            from vasp_calculation import VaspCalc
            calculation_method = VaspCalc
        elif self.create == "Gap":
            from gap_calculation import GapCalc
            calculation_method = GapCalc

        return calculation_method


    def run(self):
        path = os.path.join(self.results_directory, self.neb_directory)
        image_dirs = self.get_image_dirs(path)
        self.image_paths = self.copy_images(image_dirs)
        self.calc_method = self.calculation_method

        return



if __name__ == "__main__":
    manual = True

    if manual:
        resdir="./"
        out_dir = "start_other_neb_calc"
        create = "Vasp"

        c = ConvertNEB(resdir, out_dir, create)
        c.run()

        image_paths = c.image_paths
        method = c.calc_method

        input_directory = "input_dir"
        output_directory = "output_dir"
        potential_directory = "input_dir"

        binary = "/appl/soft/phys/vasp/6.3.0/gcc-11.2.0/bin/vasp_std"
        ncores = 128

        vasp_input_args = {"prec" : 'Accurate',
                           "xc" : 'PBE',
                           "ibrion" : -1,
                           "algo" :  'Normal',
                           "potim" :  1.0,
                           "ediff" :  1e-6,
                           "lwave" : False,
                           "lcharg" : False,
                           "ncore" : 16,
                           "kpar" : 4
                           }


        system = "CBr"

        all_calcs = []

        for path in image_paths:
            all_calcs.append( CalculationContainer(method,
                                                   CalculationData(binary = binary,
                                                                   potential_directory = potential_directory,
                                                                   input_directory     = input_directory,
                                                                   output_directory    = output_directory,
                                                                   input_args          = vasp_input_args,
                                                                   structure           = read(path, format = "xyz"),
                                                                   ncores              = ncores,
                                                                   system              = system
                                                                   )
                                                   )
                             )





        neb_args = NEB_data(images = [ci.method for ci in all_calcs],
                            n_images = len(all_calcs),
                            climb = False
                            )


        n = CalculationContainer(NEB_interface, neb_args )

        print(n)


    ############################################################
    ###---   Converting from a previous neb calculation   ---###
    ############################################################

    # prev_neb = CalculationContainer(NEB_interface, neb_args )
    # prev_neb.run()
    # os.chdir("../")



    # print(n)


    #    n.run()



    # Now get the images into other directory






# One can obtain the neb object as before from pickle

# load_from_pickle = True
# pickle_file = "neb_args.pkl"

# if load_from_pickle:
#     neb_args = pickle.load("neb_args.pkl")
