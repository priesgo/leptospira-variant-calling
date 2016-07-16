#!/usr/bin/python
import logging
import os
import sys
import scvc.config.conf_reader as conf_reader
from scvc.pipeline import Pipeline

__author__ = 'priesgo'

"""
Wrapper for running GATKs BAM realignment pipeline.
* Inputs: 
    - BAM alignments file
    - FASTA reference
* Output: 
    - Realigned BAM
"""
class BamRealignmentWrapper(Pipeline):
    
    def __init__(self):
        # Initialize parent
        super(self.__class__,self).__init__("GATKs BAM realignment wrapper.")
        # Reads configuration properties
        self.gatk_jar=conf_reader.get_property("GATK", "third-party")
        # Parses input arguments
        self.parser.add_argument('input_bam',
                            help='Input BAM alignments file')
        self.parser.add_argument('input_reference',
                            help='Input FASTA reference file')
        self.parser.add_argument('output_bam',
                            help='Output BAM')
        # It's important to avoid parsing thef first option which is parsed at hcvc
        args = self.parser.parse_args(sys.argv[2:])
        # Read input parameters
        self.input_bam = args.input_bam
        self.input_prefix = os.path.splitext(os.path.basename(self.input_bam))[0]
        self.input_reference = args.input_reference
        self.output_bam = args.output_bam
        self.output_folder = os.path.dirname(os.path.abspath(self.output_bam))
        logging.info("Input BAM alignments : %s" % self.input_bam)
        logging.info("Input FASTA reference : %s" % self.input_reference)
        logging.info("Output BAM : %s" % self.output_bam)
        
    def build_pipeline(self):
        cmd = "java -Xmx1g -jar %(gatk_jar)s -T RealignerTargetCreator -R %(input_reference)s -I %(input_bam)s \
        -o %(output_folder)s/%(input_prefix)s.intervals --filter_mismatching_base_and_quals" % { 
                                                                                                'gatk_jar':self.gatk_jar, 
                                                                                                'input_reference':self.input_reference,
                                                                                                'input_bam':self.input_bam, 
                                                                                                'output_folder':self.output_folder, 
                                                                                                'input_prefix':self.input_prefix}
        self.add_command("GATK RalignerTargetCreator", cmd, self.output_folder)
        
        cmd = "java -Xmx4g -Djava.io.tmpdir=. -jar %(gatk_jar)s -T IndelRealigner -I %(input_bam)s \
        -R %(input_reference)s -targetIntervals %(output_folder)s/%(input_prefix)s.intervals \
        -o %(output_bam)s --filter_mismatching_base_and_quals --consensusDeterminationModel USE_READS" % {
                                                                                                          'gatk_jar':self.gatk_jar, 
                                                                                                          'input_reference':self.input_reference,
                                                                                                          'input_bam':self.input_bam,
                                                                                                          'output_bam':self.output_bam,  
                                                                                                          'output_folder':self.output_folder, 
                                                                                                          'input_prefix':self.input_prefix,
                                                                                                          'output_vcf':self.output_bam}
        self.add_command("GATK IndelRealigner", cmd, self.output_folder)
        
        self.add_temporary_file("%s/%s.intervals" % (self.output_folder, self.input_prefix))
