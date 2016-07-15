#!/usr/bin/python
import logging
import os
import sys
import scvc.config.conf_reader as conf_reader
from scvc.pipeline import Pipeline

__author__ = 'priesgo'

"""
Wrapper for VEP registering a reference.
* Inputs: 
    - FASTA reference
    - GFF/GTF annotations file
    - VEP species name
* Output: 
    - The reference is registered in VEP
    
"""
class ReferenceRegisterWrapper(Pipeline):
    
    def __init__(self):
        # Initialize parent
        super(self.__class__,self).__init__("Annotation reference register wrapper.")
        # Reads configuration properties
        self.vep_home=conf_reader.get_property("VEP_HOME", "third-party")
        # Parses input arguments
        self.parser.add_argument('input_gff',
                            help='Input GFF/GTF annotations file')
        self.parser.add_argument('input_reference',
                            help='Input FASTA reference file')
        self.parser.add_argument('--species', dest='species', action='store',
                                 required=True,
                                 type=str, 
                                 help='The species reference name that will be registered in VEP')
        
        # It's important to avoid parsing thef first option which is parsed at hcvc
        args = self.parser.parse_args(sys.argv[2:])
        # Read input parameters
        self.input_gff = args.input_gff
        self.input_reference = args.input_reference
        self.species = args.species
        logging.info("Input GFF/GTF annotations: %s" % self.input_gff)
        logging.info("Input FASTA reference : %s" % self.input_reference)
        logging.info("VEP species : %s" % self.species)
        
        
    def build_pipeline(self):
        
        cmd = "perl %(vep_home)s/gtf2vep.pl -i %(input_gff)s \
        -f %(input_reference)s -d 84 -s %(species)s" % {
                                                            'vep_home':self.vep_home, 
                                                            'input_reference':self.input_reference,
                                                            'input_gff':self.input_gff,
                                                            'species':self.species}
        self.add_command("VEP reference register", cmd)