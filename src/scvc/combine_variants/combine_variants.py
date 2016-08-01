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
                                 help='Input FASTA reference file'
                                 )
        self.parser.add_argument('output_vcf',
                                 help='Output VCF file'
                                 )
        self.parser.add_argument('-V', action='append', dest='vcfs',
                                default=[], required=True,
                                help='Named VCF file to combine. Set name and file name separated by a colon (NAME:yourfile.vcf)',
                                )
        self.parser.add_argument('--method', dest='method', action='store',
                                 required=True, type=str,
                                 choices=['union', 'intersection'],
                                 help='Consensus method. One from ["union", "intersection"]'
                                 )
        
        # It's important to avoid parsing thef first option which is parsed at hcvc
        args = self.parser.parse_args(sys.argv[2:])
        # Read input parameters
        self.input_reference = args.input_reference
        self.output_vcf = args.output_vcf
        self.output_folder = os.path.dirname(os.path.abspath(self.output_vcf))
        self.output_prefix = os.path.splitext(os.path.basename(self.output_vcf))[0]
        vcfs = args.vcfs
        self.method = args.method
        self.vcfs_to_combine = {}
        for vcf in vcfs:
            values = vcf.split(":")
            if len(values) != 2:
                raise argparse.ArgumentTypeError("VCF file must be named (eg: NAME:yourfile.vcf)")
            self.vcfs_to_combine[values[0]] = values[1]
            
        logging.info("Input VCFs : %s" % str(vcfs))
        logging.info("Input FASTA reference : %s" % self.input_reference)
        logging.info("Output VCF file : %s" % self.output_vcf)
        logging.info("Consensus method : %s" % self.method)
    
      
    def build_pipeline(self):
        
        if self.method == "union":
            self.combine_variants_union()
        elif self.method == "intersection":
            self.combine_variants_intersection()
        else:
            raise NotImplemented("Consensus method for variants not implemented.")
    
    
    def combine_variants_union(self):
        
        gatk_parameters = ""
        for name, vcf in self.vcfs_to_combine.iteritems():
            gatk_parameters += " -V:%s %s" % (name, vcf)
        cmd = "java -jar %(gatk_jar)s -T CombineVariants -R %(input_reference)s %(gatk_parameters)s \
        -o %(output_vcf)s --filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED \
        --genotypemergeoption UNIQUIFY" % {
                                           'gatk_jar':self.gatk_jar, 
                                           'input_reference':self.input_reference,
                                           'output_vcf':self.output_vcf,
                                           'gatk_parameters':gatk_parameters 
                                           }
        self.add_command("GATK CombineVariants", cmd)
        
        
    def combine_variants_intersection(self):    
        
        gatk_parameters = ""
        for name, vcf in self.vcfs_to_combine.iteritems():
            gatk_parameters += " -V:%s %s" % (name, vcf)
        union_vcf = "%s/%s.tmp.union.vcf" % (self.output_folder, self.output_prefix)
        cmd = "java -jar %(gatk_jar)s -T CombineVariants -R %(input_reference)s %(gatk_parameters)s \
        -o %(union_vcf)s --filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED \
        --genotypemergeoption UNIQUIFY" % {
                                           'gatk_jar':self.gatk_jar, 
                                           'input_reference':self.input_reference,
                                           'union_vcf':union_vcf,
                                           'gatk_parameters':gatk_parameters 
                                           }
        self.add_command("GATK CombineVariants", cmd)
        
        cmd = "java -jar %(gatk_jar)s -T SelectVariants -R %(input_reference)s -V:variant %(union_vcf)s \
        -select 'set == \"Intersection\";' -o %(output_vcf)s \
        --setFilteredGtToNocall" % {
                                    'gatk_jar':self.gatk_jar, 
                                    'input_reference':self.input_reference,
                                    'union_vcf':union_vcf,
                                    'output_vcf':self.output_vcf
                                    }
        self.add_command("GATK SelectVariants", cmd)
        
        self.add_temporary_file(union_vcf + "*")  # deletes temporary union VCF
