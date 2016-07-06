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
Wrapper for running GATKs HaplotypeCaller pipeline.
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
    input_prefix = os.path.splitext(os.path.basename(input_bam))[0]
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
        
        
    logging.info("GATK PrintReads to calculate BAQ...")
    cmd = "java -jar %(gatk_jar)s -T PrintReads -R %(input_reference)s -I %(input_bam)s -baq RECALCULATE \
    -o %(output_folder)s/%(input_prefix)s.baq.bam" % {'gatk_jar':gatk_jar, 'input_reference':input_reference,
                                                      'input_bam':input_bam, 'output_folder':output_folder, 
                                                      'input_prefix':input_prefix}
    is_ok = runner.run_command(cmd=cmd, working_folder=output_folder)
    if not is_ok:
        logging.error("Error executing GATK PrintReads")
        sys.exit(1)

    # Haplotype caller variant calling pipeline
    logging.info("GATK HaplotypeCaller...")
    cmd = "java -jar %(gatk_jar)s -T HaplotypeCaller -R %(input_reference)s -I %(output_folder)s/%(input_prefix)s.baq.bam \
    --genotyping_mode DISCOVERY -stand_emit_conf 30 -stand_call_conf 30 --min_base_quality_score 13 \
    --downsampling_type NONE -ploidy %(ploidy)s -nda -allowNonUniqueKmersInRef -bamout %(output_folder)s/%(input_prefix)s.hc_reassembly.bam \
    -o %(output_vcf)s --annotateNDA --annotation BaseQualityRankSumTest --annotation ClippingRankSumTest \
    --annotation Coverage --annotation FisherStrand --annotation GCContent --annotation HomopolymerRun \
    --annotation LikelihoodRankSumTest --annotation NBaseCount --annotation QualByDepth --annotation RMSMappingQuality \
    --annotation StrandOddsRatio --annotation TandemRepeatAnnotator --annotation DepthPerAlleleBySample \
    --annotation DepthPerSampleHC --annotation StrandAlleleCountsBySample --annotation StrandBiasBySample \
    --excludeAnnotation HaplotypeScore --excludeAnnotation InbreedingCoeff" % {'gatk_jar':gatk_jar, 'input_reference':input_reference, 
                                                                               'output_folder':output_folder, 'input_prefix':input_prefix,
                                                                               'output_vcf':output_vcf, 'ploidy':ploidy}
    is_ok = runner.run_command(cmd=cmd, working_folder=output_folder)
    if not is_ok:
        logging.error("Error executing GATK HaplotypeCaller")
        sys.exit(1)
    
    # delete temporary files
    files.delete_files("%s/%s.baq.ba*" % (output_folder, input_prefix))
    
        
    logging.info("Finished GATKs HaplotypeCaller")
        
if __name__ == "__main__":
    main()
