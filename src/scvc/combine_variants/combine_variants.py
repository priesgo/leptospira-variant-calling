#!/usr/bin/python
import logging
import os
import sys
import scvc.config.conf_reader as conf_reader
from scvc.pipeline import Pipeline
import argparse

__author__ = 'priesgo'

"""
Wrapper for combining variants from several VCFs.
* Inputs: 
    - FASTA reference
    - Multiple VCF files to combine
* Output: 
    - Union VCF
    - Intersection VCF
"""
class CombineVariantsWrapper(Pipeline):
    
    def __init__(self):
        # Initialize parent
        super(self.__class__,self).__init__("Combine variants wrapper.")
        # Reads configuration properties
        self.gatk_jar=conf_reader.get_property("GATK", "third-party")
        # Parses input arguments
        self.parser.add_argument('input_reference',
                            help='Input FASTA reference file')
        self.parser.add_argument('output_union_vcf',
                            help='Output union VCF file')
        self.parser.add_argument('output_intersection_vcf',
                            help='Output intersection VCF file')
        self.parser.add_argument('-V', action='append', dest='vcfs',
                    default=[], required=True,
                    help='Named VCF file to combine. Set name and file name separated by a colon (NAME:yourfile.vcf)',
                    )
        
        # It's important to avoid parsing thef first option which is parsed at hcvc
        args = self.parser.parse_args(sys.argv[2:])
        # Read input parameters
        self.input_reference = args.input_reference
        self.output_union_vcf = args.output_union_vcf
        self.output_intersection_vcf = args.output_intersection_vcf
        vcfs = args.vcfs
        self.vcfs_to_combine = {}
        for vcf in vcfs:
            values = vcf.split(":")
            if len(values) != 2:
                raise argparse.ArgumentTypeError("VCF file must be named (eg: NAME:yourfile.vcf)")
            self.vcfs_to_combine[values[0]] = values[1]
            
        logging.info("Input VCFs : %s" % str(vcfs))
        logging.info("Input FASTA reference : %s" % self.input_reference)
        logging.info("Output union VCF file : %s" % self.output_union_vcf)
        logging.info("Output intersection VCF file : %s" % self.output_intersection_vcf)
        
    def build_pipeline(self):
        gatk_parameters = ""
        for name, vcf in self.vcfs_to_combine.iteritems():
            gatk_parameters += " -V:%s %s" % (name, vcf)
        cmd = "java -jar %(gatk_jar)s -T CombineVariants -R %(input_reference)s %(gatk_parameters)s \
        -o %(output_union_vcf)s --filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED \
        --genotypemergeoption UNIQUIFY" % {
                                           'gatk_jar':self.gatk_jar, 
                                           'input_reference':self.input_reference,
                                           'output_union_vcf':self.output_union_vcf,
                                           'gatk_parameters':gatk_parameters 
                                           }
        self.add_command("GATK CombineVariants", cmd)
        
        cmd = "java -jar %(gatk_jar)s -T SelectVariants -R %(input_reference)s -V:variant %(output_union_vcf)s \
        -select 'set == \"Intersection\";' -o %(output_intersection_vcf)s \
        --setFilteredGtToNocall" % {
                                    'gatk_jar':self.gatk_jar, 
                                    'input_reference':self.input_reference,
                                    'output_union_vcf':self.output_union_vcf,
                                    'output_intersection_vcf':self.output_intersection_vcf
                                    }
        self.add_command("GATK SelectVariants", cmd)
