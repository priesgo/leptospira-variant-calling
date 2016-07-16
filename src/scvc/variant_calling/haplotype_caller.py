#!/usr/bin/python
import logging
import os
import sys
import scvc.config.conf_reader as conf_reader
from scvc.pipeline import Pipeline

__author__ = 'priesgo'

"""
Wrapper for running GATKs HaplotypeCaller pipeline.
* Inputs: 
    - BAM alignments file
    - FASTA reference
    - Output folder
* Output: 
    - VCF
"""
class HaplotypeCallerWrapper(Pipeline):
    
    def __init__(self):
        # Initialize parent
        super(self.__class__,self).__init__("GATKs HaplotypeCaller wrapper.")
        # Reads configuration properties
        self.gatk_jar=conf_reader.get_property("GATK", "third-party")
        # Parses input arguments
        self.parser.add_argument('input_bam',
                            help='Input BAM alignments file')
        self.parser.add_argument('input_reference',
                            help='Input FASTA reference file')
        self.parser.add_argument('output_vcf',
                            help='Output VCF')
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
        self.output_folder = os.path.dirname(os.path.abspath(self.output_vcf))
        self.ploidy = str(args.ploidy)
        logging.info("Input BAM alignments : %s" % self.input_bam)
        logging.info("Input FASTA reference : %s" % self.input_reference)
        logging.info("Output VCF : %s" % self.output_vcf)
        logging.info("Ploidy : %s" % self.ploidy)
        
    def build_pipeline(self):
        cmd = "java -jar %(gatk_jar)s -T PrintReads -R %(input_reference)s -I %(input_bam)s -baq RECALCULATE \
        -o %(output_folder)s/%(input_prefix)s.baq.bam" % {'gatk_jar':self.gatk_jar, 'input_reference':self.input_reference,
                                                          'input_bam':self.input_bam, 'output_folder':self.output_folder, 
                                                          'input_prefix':self.input_prefix}
        self.add_command("GATK PrintReads to calculate BAQ", cmd, self.output_folder)
        
        cmd = "java -jar %(gatk_jar)s -T HaplotypeCaller -R %(input_reference)s -I %(output_folder)s/%(input_prefix)s.baq.bam \
        --genotyping_mode DISCOVERY -stand_emit_conf 30 -stand_call_conf 30 --min_base_quality_score 13 \
        --downsampling_type NONE -ploidy %(ploidy)s -nda -allowNonUniqueKmersInRef -bamout %(output_folder)s/%(input_prefix)s.hc_reassembly.bam \
        -o %(output_vcf)s --annotateNDA --annotation BaseQualityRankSumTest --annotation ClippingRankSumTest \
        --annotation Coverage --annotation FisherStrand --annotation GCContent --annotation HomopolymerRun \
        --annotation LikelihoodRankSumTest --annotation NBaseCount --annotation QualByDepth --annotation RMSMappingQuality \
        --annotation StrandOddsRatio --annotation TandemRepeatAnnotator --annotation DepthPerAlleleBySample \
        --annotation DepthPerSampleHC --annotation StrandAlleleCountsBySample --annotation StrandBiasBySample \
        --excludeAnnotation HaplotypeScore --excludeAnnotation InbreedingCoeff" % {'gatk_jar':self.gatk_jar, 'input_reference':self.input_reference, 
                                                                                   'output_folder':self.output_folder, 'input_prefix':self.input_prefix,
                                                                                   'output_vcf':self.output_vcf, 'ploidy':self.ploidy}
        self.add_command("GATK HaplotypeCaller", cmd, self.output_folder)
        
        self.add_temporary_file("%s/%s.baq.ba*" % (self.output_folder, self.input_prefix))
