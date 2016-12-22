#!/usr/bin/python
import logging
import os
import sys
import scvc.config.conf_reader as conf_reader
from scvc.pipeline import Pipeline

__author__ = 'priesgo'

"""
Wrapper for running Picard's BAM preprocessing pipeline.
* Inputs: 
    - BAM alignments file
* Output: 
    - Preprocessed BAM
    
#TODO: AddOrReplaceReadGroups hard codes the library and platform. parametrize it 
#TODO: parametrize duplicates removal to support PCR-free enrichment (ie: haloplex)
"""

class BamPreprocessingWrapper(Pipeline):
    
    def __init__(self):
        # Initialize parent
        super(self.__class__,self).__init__("Picard's BAM preprocessing wrapper.")
        # Reads configuration properties
        self.gatk_jar=conf_reader.get_property("PICARD", "third-party")
        self.tmp_dir=conf_reader.get_property("TMP_DIR", "third-party")
        # Parses input arguments
        self.parser.add_argument('input_bam',
                            help='Input BAM alignments file')
        self.parser.add_argument('output_bam',
                            help='Output preprocessed BAM')
        # It's important to avoid parsing thef first option which is parsed at hcvc
        args = self.parser.parse_args(sys.argv[2:])
        # Read input parameters
        self.input_bam = args.input_bam
        self.input_prefix = os.path.splitext(os.path.basename(self.input_bam))[0]
        self.output_bam = args.output_bam
        self.output_folder = os.path.dirname(os.path.abspath(self.output_bam))
        logging.info("Input BAM alignments : %s" % self.input_bam)
        logging.info("Output BAM : %s" % self.output_bam)
        
    def build_pipeline(self):
        
        cmd = "java -jar %(picard_jar)s CleanSam I=%(input_bam)s O=%(output_folder)s/%(input_prefix)s.cleaned.bam \
        VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true TMP_DIR=%(tmp_dir)s" % {'picard_jar':self.gatk_jar, 
                                                                                'input_bam':self.input_bam,
                                                                                'output_bam':self.output_bam, 
                                                                                'output_folder':self.output_folder, 
                                                                                'input_prefix':self.input_prefix,
                                                                                'tmp_dir':self.tmp_dir}
        self.add_command("Picard's CleanSam", cmd, self.output_folder)
        
        cmd = "java -jar %(picard_jar)s AddOrReplaceReadGroups I=%(output_folder)s/%(input_prefix)s.cleaned.bam \
        O=%(output_folder)s/%(input_prefix)s.readgroups.bam LB=Library PL=Illumina PU=Barcode SM=%(input_prefix)s \
        VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true TMP_DIR=%(tmp_dir)s" % {
                                                                             'picard_jar':self.gatk_jar, 
                                                                             'input_bam':self.input_bam,
                                                                             'output_bam':self.output_bam, 
                                                                             'output_folder':self.output_folder, 
                                                                             'input_prefix':self.input_prefix,
                                                                             'tmp_dir':self.tmp_dir
                                                                             }
        self.add_command("Picard's AddOrReplaceReadGroups", cmd, self.output_folder)
        
        cmd = "java -jar %(picard_jar)s SortSam I=%(output_folder)s/%(input_prefix)s.readgroups.bam \
        O=%(output_folder)s/%(input_prefix)s.readgroups.sorted.bam VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true \
        SORT_ORDER=coordinate TMP_DIR=%(tmp_dir)s" % {
                                                   'picard_jar':self.gatk_jar, 
                                                   'input_bam':self.input_bam,
                                                   'output_bam':self.output_bam, 
                                                   'output_folder':self.output_folder, 
                                                   'input_prefix':self.input_prefix,
                                                   'tmp_dir':self.tmp_dir
                                                   }
        self.add_command("Picard's SortSam", cmd, self.output_folder)
        
        cmd = "java -jar %(picard_jar)s MarkDuplicates I=%(output_folder)s/%(input_prefix)s.readgroups.sorted.bam \
        O=%(output_folder)s/%(input_prefix)s.dedupped.bam METRICS_FILE=%(output_folder)s/%(input_prefix)s.dedupped.bam.metrics \
        REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true TMP_DIR=%(tmp_dir)s \
        PROGRAM_RECORD_ID=null" % {
                                   'picard_jar':self.gatk_jar, 
                                   'input_bam':self.input_bam,
                                   'output_bam':self.output_bam, 
                                   'output_folder':self.output_folder, 
                                   'input_prefix':self.input_prefix,
                                   'tmp_dir':self.tmp_dir
                                   }
        self.add_command("Picard's MarkDuplicates", cmd, self.output_folder)
        
        cmd = "java -jar %(picard_jar)s FixMateInformation I=%(output_folder)s/%(input_prefix)s.dedupped.bam \
        O=%(output_bam)s VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true ASSUME_SORTED=false \
        TMP_DIR=%(tmp_dir)s" % {
                             'picard_jar':self.gatk_jar, 
                             'input_bam':self.input_bam,
                             'output_bam':self.output_bam, 
                             'output_folder':self.output_folder, 
                             'input_prefix':self.input_prefix,
                             'tmp_dir':self.tmp_dir
                             }
        self.add_command("Picard's FixMateInformation", cmd, self.output_folder)
        
        self.add_temporary_file("%s/%s.cleaned.ba*" % (self.output_folder, self.input_prefix))
        self.add_temporary_file("%s/%s.readgroups.ba*" % (self.output_folder, self.input_prefix))
        self.add_temporary_file("%s/%s.readgroups.sorted.ba*" % (self.output_folder, self.input_prefix))
        self.add_temporary_file("%s/%s.dedupped.ba*" % (self.output_folder, self.input_prefix))
        self.add_temporary_file("%s/snappy*.so" % (self.output_folder))
