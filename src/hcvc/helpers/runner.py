import subprocess
import shlex
import logging
import time

__author__ = 'priesgo'

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