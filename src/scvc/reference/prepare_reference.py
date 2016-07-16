#!/usr/bin/python
import logging
import os
import sys
import scvc.config.conf_reader as conf_reader
from scvc.pipeline import Pipeline

__author__ = 'priesgo'

"""
Wrapper for preparing a genomic reference for Picard and GATK.
Creates the *.dict and *.faidx from the input reference fasta.
* Inputs: 
    - FASTA reference
"""
class PrepareReferenceWrapper(Pipeline):
    
    def __init__(self):
        # Initialize parent
        super(self.__class__,self).__init__("Prepare reference genome wrapper.")
        # Reads configuration properties
        self.samtools_home=conf_reader.get_property("SAMTOOLS_HOME", "third-party")
        self.picard_jar=conf_reader.get_property("PICARD", "third-party")
        # Parses input arguments
        self.parser.add_argument('input_reference',
                            help='Input FASTA reference file')
        # It's important to avoid parsing thef first option which is parsed at hcvc
        args = self.parser.parse_args(sys.argv[2:])
        # Read input parameters
        self.input_reference = args.input_reference
        self.input_prefix = os.path.splitext(os.path.basename(self.input_reference))[0]
        self.output_folder = os.path.dirname(os.path.abspath(self.input_reference))                                             
        logging.info("Input FASTA reference : %s" % self.input_reference)
        
    def build_pipeline(self):
        cmd = "%(samtools_home)s/samtools faidx %(input_reference)s" % {'samtools_home':self.samtools_home, 
                                                                        'input_reference':self.input_reference}
        self.add_command("Index fasta file", cmd)
        
        cmd = "java -jar %(picard_jar)s CreateSequenceDictionary R=%(input_reference)s \
        O=%(output_folder)s/%(input_prefix)s.dict" % {'picard_jar':self.picard_jar, 
                                                      'input_reference':self.input_reference, 
                                                      'output_folder':self.output_folder, 
                                                      'input_prefix':self.input_prefix}
        self.add_command("Create Picard reference dictionary", cmd, self.output_folder)
