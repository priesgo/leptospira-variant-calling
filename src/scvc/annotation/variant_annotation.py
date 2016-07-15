#!/usr/bin/python
import logging
import os
import sys
import scvc.config.conf_reader as conf_reader
from scvc.pipeline import Pipeline

__author__ = 'priesgo'

"""
Wrapper for VEP variant annotation.
* Inputs: 
    - VCF variants file
    - FASTA reference
    - VEP species name
* Output: 
    - Annotated VCF
    
"""
class VariantAnnotationWrapper(Pipeline):
    
    def __init__(self):
        # Initialize parent
        super(self.__class__,self).__init__("Variant annotation wrapper.")
        # Reads configuration properties
        self.vep_home=conf_reader.get_property("VEP_HOME", "third-party")
        # Parses input arguments
        self.parser.add_argument('input_vcf',
                            help='Input VCF variants file')
        self.parser.add_argument('input_reference',
                            help='Input FASTA reference file')
        self.parser.add_argument('output_vcf',
                            help='Output filtered VCF')
        self.parser.add_argument('--species', dest='species', action='store',
                                 required=True,
                                 type=str, 
                                 help='The species reference name as registered in VEP')
        
        # It's important to avoid parsing thef first option which is parsed at hcvc
        args = self.parser.parse_args(sys.argv[2:])
        # Read input parameters
        self.input_vcf = args.input_vcf
        self.input_prefix = os.path.splitext(os.path.basename(self.input_vcf))[0]
        self.input_reference = args.input_reference
        self.output_vcf = args.output_vcf
        self.output_folder = os.path.dirname(self.output_vcf)
        self.species = args.species
        logging.info("Input VCF alignments : %s" % self.input_vcf)
        logging.info("Input FASTA reference : %s" % self.input_reference)
        logging.info("Output VCF : %s" % self.output_vcf)
        logging.info("VEP species : %s" % self.species)
        
        
    def build_pipeline(self):
        
        cmd = "perl %(vep_home)s/variant_effect_predictor.pl -i %(input_vcf)s \
        -o %(output_vcf)s --species %(species)s \
        --vcf --offline --symbol --protein --uniprot --biotype --hgvs \
        --fasta %(input_reference)s --domains --total_length \
        --variant_class --pick --no_escape" % {
                                                            'vep_home':self.vep_home, 
                                                            'input_reference':self.input_reference,
                                                            'input_vcf':self.input_vcf,
                                                            'output_vcf':self.output_vcf,
                                                            'species':self.species}
        self.add_command("VEP variant annotation", cmd, self.output_folder)