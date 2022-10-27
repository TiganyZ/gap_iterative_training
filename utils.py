#!/usr/bin/python3
import os, subprocess, shutil
from glob import glob


class Utils:

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




    
    def wrap_function(self, prefix, function, message):
        print(f"> {prefix}:          {message} ")
        function()
        print(f"> {prefix}: SUCCESS: {message} \n")

        return True


    def check_key(self, dic, key):
        ret = False
        if key in dic.keys():
            ret = True
        else:
            ret = False
        return ret



    def check_keys(self, dic, keys = ( "input_directory", "potential_directory", "output_directory" )):

        for key in keys:
            if not self.check_key(dic, key):
                print(f"""\n{dic}\n\n
                ########################################################
                ###---   WARNING! The key = {key} is not found!   ---###
                ########################################################\n   --> Exiting <--\n""")
                exit(1)



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
            print(f"""
            -->   Found file {file} in ./ directory    <--
            """)
        elif self.check_file(f"{dir}/{subdir}", file):
            print(f"""
            -->   Found file {file} in {subdir} subdirectory    <--
            """)
            file = f"{dir}/{subdir}/{file}"
        else:
            print(f"""
                ###############################################################################
                ###---   FATAL: No file {file} in directory or {subdir} subdirectory.    ---###
                ###############################################################################
                """)
            return None

        return file



    def check_file(self, directory, file, stop=False):
        if not os.path.exists(f"{directory}/{file}"):
            print(f"""
            ################################################################
            ###---   WARNING! There is no path {directory}/{file} !   ---###
            ################################################################\n""")
            if stop:
                print("\n  --> Rectify this <--\n  I'm leaving you...")
                exit(1)
            else:
                return False
        else:
            return True

    def check_copy_file(self, directory, file, output_directory):
        if self.check_file(directory, file) and os.path.exists(f"{output_directory}"):
            shutil.copy( f"{directory}/{file}", f"{output_directory}/" )

    def save_file_in_dir(self, file, dir, savedir ):
        if os.path.exists(f"{dir}/{file}"):
            if not os.path.exists(f"{dir}/{savedir}"):
                os.mkdir(f"{dir}/{savedir}")
            n_files = len(os.listdir(f"{dir}/{savedir}"))
            if n_files > 0:
                # Check if the file is the same or not
                cmd = "diff {dir}/{file} {dir}/{savedir}/{file}_{n_files - 1} | wc -l"
                l = self.utils.piped_subprocess(cmd)
                if int(l) > 0:
                    f = file.split(".")
                    if len(f) > 1:
                        f[-1] = f"_{n_files}." + f[-1]
                        filename = ''.join(f)
                    else:
                        filename = f"{file}_{n_files}"
                    shutil.copy(f"{dir}/{file}", f"{dir}/{savedir}/{filename}")


    def piped_subprocess(self, commands, file=None):
        for i, cmd in enumerate(commands.split("|")):
            if i == 0:
                p = subprocess.Popen(cmd.split(), stdout = subprocess.PIPE )
            if i == len(commands.split("|"))-1 and file is not None:
                p = subprocess.Popen(cmd.split(), stdin=p.stdout, stdout=file  )
            else:
                p = subprocess.Popen(cmd.split(), stdin=p.stdout, stdout = subprocess.PIPE )

            self.check_subprocess(p)

        out, errs = p.communicate()
        return out

    def check_subprocess(self, out, verbosity=100):

        if out.returncode != 0:
            print(f"""
            ##################################################################################
            ###---   WARNING! Something has gone wrong with the subprocess {out.args}   ---###
            ##################################################################################
            """)
            exit(1)

        elif verbosity > 50:
            print(f"CHECK_SUBPROCESS: {out}")


    def get_save_name(self, path, result, prefix, ext=".json"):
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
