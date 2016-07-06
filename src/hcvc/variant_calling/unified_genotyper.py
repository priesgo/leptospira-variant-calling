#!/usr/bin/python
import argparse
import logging
import os
import sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/../..")
import hcvc.helpers.runner as runner
import hcvc.helpers.files as files
import hcvc.config.conf_reader as conf_reader

__author__ = 'priesgo'

"""
Wrapper for running GATKs UnifiedGenpotyper pipeline.
* Inputs: 
    - BAM alignments file
    - FASTA reference
    - Output folder
* Output: 
    - VCF
"""
    

def main():

    # read this value from a config file
    gatk_jar=conf_reader.get_property("GATK", "third-party")
    
    logging.basicConfig(level=logging.INFO)
    
    # Configure input parameters
    parser = argparse.ArgumentParser(description="Runs GATK's Haplotype Caller.")
    parser.add_argument('input_bam',
                        help='Input BAM alignments file')
    parser.add_argument('input_reference',
                        help='Input FASTA reference file')
    parser.add_argument('output_vcf',
                        help='Output folder')
    parser.add_argument('--ploidy', dest='ploidy', action='store',
                        default=1, type=int, 
                        help='The ploidy of the organism.')
    args = parser.parse_args()
    
    # Read input parameters
    input_bam = args.input_bam
    input_reference = args.input_reference
    output_vcf = args.output_vcf
    output_folder = os.path.dirname(output_vcf)
    ploidy = str(args.ploidy)
    logging.info("Input BAM alignments : %s" % input_bam)
    logging.info("Input FASTA reference : %s" % input_reference)
    logging.info("Output VCF : %s" % output_vcf)
    logging.info("Ploidy : %s" % ploidy)
    
    # creates output folder in case it does not exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        
    logging.info("GATK UnifiedGenotyper ...")
    cmd = "java -jar %(gatk_jar)s -T UnifiedGenotyper -R %(input_reference)s -I %(input_bam)s --genotyping_mode DISCOVERY \
    -stand_emit_conf 30 -stand_call_conf 30 --min_base_quality_score 13 -baq RECALCULATE --downsampling_type NONE \
    -ploidy %(ploidy)s -nda -o %(output_vcf)s --annotateNDA --annotation BaseQualityRankSumTest \
    --annotation ClippingRankSumTest \
    --annotation Coverage --annotation FisherStrand --annotation GCContent --annotation HomopolymerRun \
    --annotation LikelihoodRankSumTest --annotation NBaseCount --annotation QualByDepth --annotation RMSMappingQuality \
    --annotation StrandOddsRatio --annotation TandemRepeatAnnotator --annotation DepthPerAlleleBySample \
    --annotation DepthPerSampleHC --annotation StrandAlleleCountsBySample --annotation StrandBiasBySample \
    --excludeAnnotation HaplotypeScore --excludeAnnotation InbreedingCoeff" % {'gatk_jar':gatk_jar, 
                                                                               'input_reference':input_reference, 
                                                                               'output_vcf':output_vcf, 
                                                                               'ploidy':ploidy}
    is_ok = runner.run_command(cmd=cmd, working_folder=output_folder)
    if not is_ok:
        logging.error("Error executing GATK UnifiedGenotyper")
        sys.exit(1)
        
    logging.info("Finished GATKs UnifiedGenotyper")
        
if __name__ == "__main__":
    main()
