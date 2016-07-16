#!/usr/bin/python
import logging
import os
import sys
import scvc.config.conf_reader as conf_reader
from scvc.pipeline import Pipeline

__author__ = 'priesgo'

"""
Wrapper for running Samtools pileup pipeline.
* Inputs: 
    - BAM alignments file
    - FASTA reference
    - Output folder
* Output: 
    - VCF
"""
class SamtoolsPileupWrapper(Pipeline):
    
    def __init__(self):
        # Initialize parent
        super(self.__class__,self).__init__("Samtools pileup wrapper.")
        # Reads configuration properties
        self.gatk_jar=conf_reader.get_property("GATK", "third-party")
        self.samtools_home=conf_reader.get_property("SAMTOOLS_HOME", "third-party")
        self.bcftools_home=conf_reader.get_property("BCFTOOLS_HOME", "third-party")
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
        self.output_folder = os.path.dirname(os.path.abspath(self.output_vcf))
        self.ploidy = str(args.ploidy)
        logging.info("Input BAM alignments : %s" % self.input_bam)
        logging.info("Input FASTA reference : %s" % self.input_reference)
        logging.info("Output VCF : %s" % self.output_vcf)
        logging.info("Ploidy : %s" % self.ploidy)
        
    def build_pipeline(self):
        cmd = "%(samtools_home)s/samtools mpileup --min-BQ 13 --adjust-MQ 50 --redo-BAQ --min-MQ 1 --illumina1.3+ \
        --output-BP --output-MQ --uncompressed --fasta-ref %(input_reference)s %(input_bam)s \
        -o %(output_folder)s/%(input_prefix)s.tmp.bcf" % {
                                                          'samtools_home':self.samtools_home,
                                                          'input_reference':self.input_reference,
                                                          'input_bam':self.input_bam, 
                                                          'output_folder':self.output_folder, 
                                                          'input_prefix':self.input_prefix}
        self.add_command("Samtools mpileup", cmd, self.output_folder)
        
        cmd = "%(bcftools_home)s/bcftools call --multiallelic-caller --variants-only --output-type v \
        --ploidy %(ploidy)s %(output_folder)s/%(input_prefix)s.tmp.bcf" % {
                                                                             'bcftools_home':self.bcftools_home,
                                                                             'output_folder':self.output_folder, 
                                                                             'input_prefix':self.input_prefix,
                                                                             'ploidy':self.ploidy}
        self.add_command("bcftools call", cmd, self.output_folder, output_file="%(output_folder)s/%(input_prefix)s.tmp.vcf" % {'output_folder':self.output_folder, 
                                                                                                                               'input_prefix':self.input_prefix})
        
        cmd = "java -jar %(gatk_jar)s -T VariantAnnotator -I %(input_bam)s -R %(input_reference)s \
        -V %(output_folder)s/%(input_prefix)s.tmp.vcf \
        -o %(output_vcf)s --annotation BaseQualityRankSumTest --annotation ClippingRankSumTest --annotation Coverage \
        --annotation FisherStrand --annotation GCContent --annotation HomopolymerRun --annotation LikelihoodRankSumTest \
        --annotation NBaseCount --annotation QualByDepth --annotation RMSMappingQuality --annotation StrandOddsRatio \
        --annotation TandemRepeatAnnotator --annotation DepthPerAlleleBySample --annotation DepthPerSampleHC \
        --annotation StrandAlleleCountsBySample --annotation StrandBiasBySample --excludeAnnotation HaplotypeScore \
        --excludeAnnotation InbreedingCoeff" % {'gatk_jar':self.gatk_jar,
                                                'input_bam':self.input_bam,
                                                'input_reference':self.input_reference, 
                                                'output_folder':self.output_folder, 
                                                'input_prefix':self.input_prefix,
                                                'output_vcf':self.output_vcf}
        self.add_command("GATK VariantAnnotator", cmd, self.output_folder)
        
        
        self.add_temporary_file("%s/%s.tmp.vcf" % (self.output_folder, self.input_prefix))
        self.add_temporary_file("%s/%s.tmp.bcf" % (self.output_folder, self.input_prefix))