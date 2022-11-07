#!/usr/bin/python3
import os, subprocess, shutil
from glob import glob
from functools import partial 


class Utils:

    def print_statement(self, message, warning = "WARNING: ", buff_char = "#"):
        # Add other characters to the message
        if warning is None:
            warning = ""
        message = buff_char*3 + "-->   " + warning + message + "   <--" + buff_char*3
        buff = buff_char * len(message)

        print("\n" + buff + "\n" + message + "\n" + buff + "\n")

    def warning(self, message):
        return self.print_statement(message, buff_char="!")
    
    def notice(self, message):
        return self.print_statement(message, warning="NOTICE: ", buff_char="~")
    
    def success(self, message):
        return self.print_statement(message, warning="SUCCESS: ", buff_char="-")

    def fatal(self, message):
        return self.print_statement(message, warning="FATAL: ", buff_char="#")
        

    def wrap_function(self, prefix, function, message):
        print(f"> {prefix}:          {message} ")
        function()
        print(f"> {prefix}: SUCCESS: {message} \n")

        return True


    ##############################
    ###---   Dictionaries   ---###
    ##############################

    def check_key(self, dic, key):
        ret = False
        if key in dic.keys():
            ret = True
        else:
            ret = False
        return ret



    def check_keys(self, dic, keys = ( "input_directory", "potential_directory", "output_directory" ), verbosity=0):
        for key in keys:
            if not self.check_key(dic, key):
                if verbosity > 50:
                    print(f"\n{dic}")
                    self.print_statement(message, f"The key = {key} is not found!")
                raise ValueError



    #####################################################
    ###---   Copying and checking for files/dirs   ---###
    #####################################################

    def check_copy_tree(self, src, dst):
        # Check that directory dst does not exist
        if not os.path.exists(dst):
            shutil.copytree(src, dst)
        else:
            ldst = glob(f"{dst}/*/", recursive = True)
            lsrc = glob(f"{src}/*/", recursive = True)

            for dir in lsrc:
                if dir in ldst:
                    print(f"check_copy_tree: {src}/{dir} has same name as {dst}/dir. Copying just the files")
                    for f in os.listdir(src):
                        if os.path.isfile(f):
                            end_dir = dir.split("/")[-1]
                            print(f"check_copy_tree: copying {src}/{f} to {dst}/{end_dir}/")
                            shutil.copy(f, f"{dst}/{end_dir}/")

    def copy_only_files(self, src, dst):
        files = os.listdir(src)
        for f in files:
            if os.path.isfile(f):
                print(f"Utils.copy_only_files: copying {f}, {dst}")
                shutil.copy(f, f"{dst}/")


    def check_file_dir_subdir(self, file, dir ='.', subdir="gap_files"):
        if self.check_file(dir, file):
            self.print_statement(f"Found file {file} in ./ directory", warning=None)

        elif self.check_file(f"{dir}/{subdir}", file):
            self.print_statement(f"Found file {file} in {subdir} subdirectory", warning=None)
            file = f"{dir}/{subdir}/{file}"
        else:
            self.print_statement(f"No file {file} in directory or {subdir} subdirectory.")

        return file



    def check_file(self, directory, file, stop=False):
        if not os.path.exists( os.path.normpath( os.path.join(directory,file) ) ):
            self.print_statement(f"There is no path {os.path.join(directory,file)}!")
            if stop:
                print("\n  --> Rectify this <--\n  I'm leaving you...")
                raise ValueError
            else:
                return False
        else:
            return True

    def check_copy_file(self, directory, file, output_directory):
        if self.check_file(directory, file) and os.path.exists(f"{output_directory}"):
            shutil.copy( os.path.join(directory,file), f"{output_directory}/" )


    #################################
    ###---   Subprocess & os   ---###
    #################################
    def piped_subprocess(self, commands, file=None):
        if file is not None:
            with open(file, "w") as fil:
                out=subprocess.run(commands,  stdout = fil, shell=True)
                self.check_subprocess(out)
                return True
        else:
            p = subprocess.Popen(commands,  stdout = subprocess.PIPE, shell=True)
            out, err = p.communicate()
            return out.decode("utf-8")
        # p = {}
        # for i, cmd in enumerate(commands.split("|")):
        #     print("piped_subprocess > ", i, cmd.split())
        #     if i == 0:
        #         p[i] = subprocess.Popen(cmd.split(), stdout = subprocess.PIPE, stderr=subprocess.PIPE )
        #     elif i == len(commands.split("|"))-1 and file is not None:
        #         p[i] = subprocess.Popen(cmd.split(), stdin=p[i-1].stdout, stdout=file  )
        #     else:
        #         p[i] = subprocess.Popen(cmd.split(), stdin=p[i-1].stdout, stdout = subprocess.PIPE )

        # for i in reversed(range(len(commands.split("|")))):
        #     out, errs = p[i].communicate()

        # self.check_subprocess(p[0], verbosity=100)
        # return out


    def check_subprocess(self, out, verbosity=50):
        print(out)

        if out.returncode != 0:
            self.fatal(f"Something has gone wrong with the subprocess \"{out.args}\"")
            raise RuntimeError

        elif verbosity > 50:
            print(f"CHECK_SUBPROCESS: {out}")


    ########################
    ###---   Saving   ---###
    ########################

    def get_save_name(self, path, result, prefix, ext=".json"):
        if not os.path.exists(path):
            os.mkdir(path)

        l = os.listdir(path)

        if len(prefix) == 0:
            name = '_'.join( list(result.keys()) )
        else:
            name = prefix
        suffix=''
        counter = 0
        while name+suffix+ext in l:
            suffix = f"_{counter:d}"
            counter += 1

        return name+suffix+ext


    ##################################################################################
    ###---   Generating driver file if necessary for some future calculations   ---###
    ##################################################################################

    def get_driver_template(self, driver_args, just_command=True, job_blurb=None):

        if self.check_key(driver_args, "modules"):
            driver_args["all_modules"] = "\n ".join( [ "module load " + mod for mod in driver_args["modules"] ] )
        else:
            driver_args["all_modules"] = ""

        if self.check_key(driver_args, "export_paths"):
            driver_args["all_exports"] = "\n ".join( [ "export " + mod      for mod in driver_args["export_paths"] ] )
        else:
            driver_args["all_exports"] = ""

        if job_blurb is None:
            driver_template = """
            #!/bin/bash -l
            #SBATCH --job-name=%(job-name)s
            #SBATCH --account=%(account)s
            #SBATCH -p %(queue)s
            #SBATCH --ntasks-per-node=%(ntasks-per-node)s
            #SBATCH --nodes=%(nodes)s
            #SBATCH -t %(walltime)s
            #SBATCH -o %(output)s

            module reset
            %(all_modules)s
            set -xe

            %(all_exports)s

            %(command)s -n %(ncores)s %(binary)s

            """ % ( driver_args )
        else:
            driver_template = job_blurb + """

            module reset
            %(all_modules)s
            set -xe

            %(all_exports)s

            %(command)s -n %(ncores)s %(binary)s

            """
            driver_template = driver_template % ( driver_args )


        if just_command:
            driver_template = """
            #!/bin/bash
            %(all_modules)s

            %(all_exports)s

            %(command)s -n %(ncores)s %(binary)s

            """ % ( driver_args )


        return driver_template
