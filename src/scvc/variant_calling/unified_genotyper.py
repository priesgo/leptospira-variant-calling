#!/usr/bin/python
import logging
import os
import sys
import scvc.config.conf_reader as conf_reader
from scvc.pipeline import Pipeline

__author__ = 'priesgo'

"""
Wrapper for running GATKs UnifiedGenotyper pipeline.
* Inputs: 
    - BAM alignments file
    - FASTA reference
    - Output folder
* Output: 
    - VCF
"""
class UnifiedGenotyperWrapper(Pipeline):
    
    def __init__(self):
        # Initialize parent
        super(self.__class__,self).__init__("GATKs UnifiedGenotyper wrapper.")
        # Reads configuration properties
        self.gatk_jar=conf_reader.get_property("GATK", "third-party")
        # Parses input arguments
        self.parser.add_argument('input_bam',
                            help='Input BAM alignments file')
        self.parser.add_argument('input_reference',
                            help='Input FASTA reference file')
        self.parser.add_argument('output_vcf',
                            help='Output folder')
        self.parser.add_argument('--ploidy', dest='ploidy', action='store',
                            default=1, type=int, 
                            help='The ploidy of the organism.')
        # It's important to avoid parsing thef first option which is parsed at hcvc
        args = self.parser.parse_args(sys.argv[2:])
        # Read input parameters
        self.input_bam = args.input_bam
        self.input_prefix = os.path.splitext(os.path.basename(self.input_bam))[0]
        self.input_reference = args.input_reference
        self.output_vcf = args.output_vcf
        self.output_folder = os.path.dirname(self.output_vcf)
        self.ploidy = str(args.ploidy)
        logging.info("Input BAM alignments : %s" % self.input_bam)
        logging.info("Input FASTA reference : %s" % self.input_reference)
        logging.info("Output VCF : %s" % self.output_vcf)
        logging.info("Ploidy : %s" % self.ploidy)
        
    def build_pipeline(self):
        cmd = "java -jar %(gatk_jar)s -T UnifiedGenotyper -R %(input_reference)s -I %(input_bam)s --genotyping_mode DISCOVERY \
        -stand_emit_conf 30 -stand_call_conf 30 --min_base_quality_score 13 -baq RECALCULATE --downsampling_type NONE \
        -ploidy %(ploidy)s -nda -o %(output_vcf)s --annotateNDA --annotation BaseQualityRankSumTest \
        --annotation ClippingRankSumTest \
        --annotation Coverage --annotation FisherStrand --annotation GCContent --annotation HomopolymerRun \
        --annotation LikelihoodRankSumTest --annotation NBaseCount --annotation QualByDepth --annotation RMSMappingQuality \
        --annotation StrandOddsRatio --annotation TandemRepeatAnnotator --annotation DepthPerAlleleBySample \
        --annotation DepthPerSampleHC --annotation StrandAlleleCountsBySample --annotation StrandBiasBySample \
        --excludeAnnotation HaplotypeScore --excludeAnnotation InbreedingCoeff" % {'gatk_jar':self.gatk_jar, 
                                                                                   'input_reference':self.input_reference, 
                                                                                   'output_vcf':self.output_vcf, 
                                                                                   'ploidy':self.ploidy}
        self.add_command("GATK UnifiedGenotyper", cmd, self.output_folder)
