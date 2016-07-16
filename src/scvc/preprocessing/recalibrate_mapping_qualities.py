#!/usr/bin/python
import logging
import os
import sys
import scvc.config.conf_reader as conf_reader
from scvc.pipeline import Pipeline

__author__ = 'priesgo'

"""
Wrapper for running GATKs mapping qualities recalibration pipeline.
* Inputs: 
    - BAM alignments file
    - FASTA reference
* Output: 
    - Recalibrated BAM
"""
class RecalibrateMappingQualitiesWrapper(Pipeline):
    
    def __init__(self):
        # Initialize parent
        super(self.__class__,self).__init__("GATKs mapping qualities recalibration wrapper.")
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
        """
        java -Xmx1g -jar $GATK -T PrintReads -I $INPUT_BAM -R $REFERENCE -BQSR $OUTPUT_DIR/$PREFIX_LOCAL.recal_data.grp --out $OUTPUT_BAM --filter_mismatching_base_and_quals --allow_potentially_misencoded_quality_scores
        java -jar $GATK -T AnalyzeCovariates -R $REFERENCE -BQSR $OUTPUT_DIR/$PREFIX_LOCAL.recal_data.grp -plots $OUTPUT_DIR/$PREFIX_LOCAL.BQSR.pdf
        rm -f $OUTPUT_DIR/$PREFIX_LOCAL.recal_data.grp
        """
        
        
        cmd = "java -Xmx1g -jar %(gatk_jar)s -T PrintReads -I %(input_bam)s -R %(input_reference)s \
        -BQSR %(output_folder)s/%(input_prefix)s.recal_data.grp --out %(output_bam)s \
        --filter_mismatching_base_and_quals \
        --allow_potentially_misencoded_quality_scores" % {
                                                          'gatk_jar':self.gatk_jar, 
                                                          'input_reference':self.input_reference,
                                                          'input_bam':self.input_bam,
                                                          'output_bam':self.output_bam, 
                                                          'output_folder':self.output_folder, 
                                                          'input_prefix':self.input_prefix
                                                          }
        self.add_command("GATK PrintReads", cmd, self.output_folder)
        
        cmd = "java -jar %(gatk_jar)s -T AnalyzeCovariates -R %(input_reference)s \
        -BQSR %(output_folder)s/%(input_prefix)s.recal_data.grp \
        -plots %(output_folder)s/%(input_prefix)s.BQSR.pdf" % {
                                                      'gatk_jar':self.gatk_jar, 
                                                      'input_reference':self.input_reference,
                                                      'output_folder':self.output_folder, 
                                                      'input_prefix':self.input_prefix
                                                      }
        self.add_command("GATK AnalyzeCovariates", cmd, self.output_folder)
        
        self.add_temporary_file("%s/%s.recal_data.grp" % (self.output_folder, self.input_prefix))
