import glob
import os
import logging


def delete_files(file_name_pattern):
    for file_name in glob.glob(file_name_pattern):
        try:
            os.remove(file_name)
            logging.error("File [%s] is removed!" % (file_name))
        except Exception, e:
            logging.error("Failed to remove file [%s]: %s" % (file_name, str(e)))