import time
import shlex
import subprocess
import argparse
import logging
import sys
import os
import glob

class Pipeline(object):
    
    
    def __init__(self, description):
        # Initialize argparser
        self.parser = argparse.ArgumentParser(description=description)
        # Initialize commands
        self.commands = []
        # Initialize temp files
        self.temp_files = []
        
    def build_pipeline(self):
        """
        Abstract method
        """
        raise NotImplementedError("Please Implement this method")
        
    def add_command(self, name, cmd, output_folder=None, output_file=None):
        """
        Adds a command. Order matters.
        """
        self.commands.append({"name":name,
                              "command":cmd,
                              "output_folder":output_folder,
                              "output_file":output_file})
    
    def run_sequential_pipeline(self):
        """
        Runs a sequential pipeline.
        PRE: self.command has been preloaded
        """
        # Calls an abstract function
        self.build_pipeline()
        for command in self.commands:
            # Reads command data
            name = command["name"] if "name" in command else None
            cmd = command["command"] if "command" in command else None
            output_folder = command["output_folder"] if "output_folder" in command else None
            output_file = command["output_file"] if "output_file" in command else None
            # Creates output folder if necessary
            if output_folder is not None and not os.path.exists(output_folder):
                os.makedirs(output_folder)
            # Runs the command
            logging.info("Running %s ..." % name)
            is_ok = run_command(cmd=cmd, working_folder=output_folder, output_file=output_file)
            if not is_ok:
                logging.error("Error executing %s" % name)
                sys.exit(1)
    
    def add_temporary_file(self, temp_file):
        self.temp_files.append(temp_file)
    
    def delete_temporary_files(self):
        # delete temporary files
        for file_name_pattern in self.temp_files:
            for file_name in glob.glob(file_name_pattern):
                try:
                    os.remove(file_name)
                    logging.error("File [%s] is removed!" % (file_name))
                except Exception, e:
                    logging.error("Failed to remove file [%s]: %s" % (file_name, str(e)))
            
def run_command(cmd, working_folder=None, output_file=None):
    """
    Runs a command line in a child process and waits for it to finish. 
    Returns a boolean with return status.
    """    
    logging.info("Running command [%s]" % cmd)
    start = time.time()
    args = shlex.split(cmd)
    is_correct_execution = True
    try:
        p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=-1, cwd=working_folder)
        out, err = p.communicate()
    except Exception, e:
        is_correct_execution = False
        logging.error("Unknown error: [%s]" % str(e))
        return is_correct_execution
    
    # checks the return code
    if p.returncode != 0:
        is_correct_execution = False
        logging.error("Unknown error: [%s]" % p.returncode)
    # prints stdout and stderr if contain any info
    if out:
        logging.info("stdout: %s" % out)
        print out
    if output_file is not None:
        f = open(output_file, "w")
        f.write(out)
        f.close()
    if err:
        logging.error("stderr: %s" % err)
    end = time.time()
    logging.info("Finished running command. Elapsed time %s seconds" % str(end-start))
    
    return is_correct_execution