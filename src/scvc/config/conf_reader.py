__author__ = 'guanabana'
import os
import sys
import codecs
from ConfigParser import SafeConfigParser

class ConfReaderException(Exception):
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)


def get_config_filepath():
    config_folder = os.path.normpath(os.path.dirname(os.path.realpath(__file__)))
    config_file = "%s/scvc.ini" % (config_folder)
    return config_file


def get_property(option, section):
    config_file = get_config_filepath()
    if os.path.exists(config_file):
        parser = SafeConfigParser()
        #parser.readfp(codecs.open(config_file, "r", "utf-8"))
        parser.read(config_file)
        if not parser.has_section(section):
            raise ConfReaderException("Section '%s' does not exist" % (section))
        elif not parser.has_option(section, option):
            raise ConfReaderException("Option '%s' at section '%s' does not exist" % (option, section))
        else:
            # Gets the value as a string when it exists
            value = parser.get(section, option, raw=True)
    else:
        # Module does not exist
        raise ConfReaderException("Configuration file '%s' does not exist" % config_file)

    return value

