#!/usr/bin/python3
import os, subprocess, shutil
from glob import glob


class Utils:

    def check_key(self, dic, key):
        print(f" > checkkey: is key = \"{key}\" present?")
        ret = False
        if key in dic.keys():
            print("Present, ", end =" ")
            print("value =", dic[key])
            ret = True
        else:
            print(f" {key} Not present")
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
                    shutil.copy(f, f"{dst}/")



    def check_file(self, directory, file):
        if not os.path.exists(f"{directory}/{file}"):
            print(f"""
            ################################################################
            ###---   WARNING! There is no path {directory}/{file} !   ---###
            ################################################################\n\n  --> Rectify this <--\n  I'm leaving you...""")
            exit(1)
        else:
            return True

    def check_copy_file(self, directory, file, output_directory):
        if self.check_file(directory, file) and os.path.exists(f"{output_directory}"):
            shutil.copy( f"{directory}/{file}", f"{output_directory}/" )

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

    def check_subprocess(self, out):

        if out.returncode != 0:
            print(f"""
            ##################################################################################
            ###---   WARNING! Something has gone wrong with the subprocess {out.args}   ---###
            ##################################################################################
            """)
            exit(1)

        if self.verbosity > 50:
            print(f"CHECK_SUBPROCESS: {out}")
