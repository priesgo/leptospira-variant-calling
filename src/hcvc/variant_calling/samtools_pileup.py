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
    samtools_home=conf_reader.get_property("SAMTOOLS_HOME", "third-party")
    bcftools_home=conf_reader.get_property("BCFTOOLS_HOME", "third-party")
    
    logging.basicConfig(level=logging.INFO)
    
    # Configure input parameters
    parser = argparse.ArgumentParser(description="Runs Samtools pileup variant calling.")
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
         
    logging.info("Samtools mpileup | bcftools call...")
    cmd = "%(samtools_home)s/samtools mpileup --min-BQ 13 --adjust-MQ 50 --redo-BAQ --min-MQ 1 --illumina1.3+ \
    --output-BP --output-MQ --uncompressed --fasta-ref %(input_reference)s %(input_bam)s | \
    %(bcftools_home)s/bcftools call --multiallelic-caller --variants-only --output-type v \
    --ploidy %(ploidy)s > %(output_folder)s/%(input_prefix)s.tmp.vcf" % {'samtools_home':samtools_home,
                                                                         'bcftools_home':bcftools_home,
                                                                         'input_reference':input_reference,
                                                                         'input_bam':input_bam, 
                                                                         'output_folder':output_folder, 
                                                                         'input_prefix':input_prefix,
                                                                         'ploidy':ploidy}
    is_ok = runner.run_command(cmd=cmd, working_folder=output_folder)
    if not is_ok:
        logging.error("Error executing Samtools pileup")
        sys.exit(1)

    # Haplotype caller variant calling pipeline
    logging.info("GATK VariantAnnotator...")
    cmd = "java -jar %(gatk_jar)s -T VariantAnnotator -I %(input_bam)s -R %(input_reference)s \
    -V %(output_folder)s/%(input_prefix)s.tmp.vcf \
    -o %(output_vcf)s --annotation BaseQualityRankSumTest --annotation ClippingRankSumTest --annotation Coverage \
    --annotation FisherStrand --annotation GCContent --annotation HomopolymerRun --annotation LikelihoodRankSumTest \
    --annotation NBaseCount --annotation QualByDepth --annotation RMSMappingQuality --annotation StrandOddsRatio \
    --annotation TandemRepeatAnnotator --annotation DepthPerAlleleBySample --annotation DepthPerSampleHC \
    --annotation StrandAlleleCountsBySample --annotation StrandBiasBySample --excludeAnnotation HaplotypeScore \
    --excludeAnnotation InbreedingCoeff" % {'gatk_jar':gatk_jar,
                                            'input_bam':input_bam,
                                            'input_reference':input_reference, 
                                            'output_folder':output_folder, 
                                            'input_prefix':input_prefix,
                                            'output_vcf':output_vcf}
    is_ok = runner.run_command(cmd=cmd, working_folder=output_folder)
    if not is_ok:
        logging.error("Error executing GATK VariantAnnotator")
        sys.exit(1)
    
    # delete temporary files
    files.delete_files("%(output_folder)s/%(input_prefix)s.tmp.vcf" % {'output_folder':output_folder, 
                                                                       'input_prefix':input_prefix})
    
        
    logging.info("Finished Samtools variant calling")
        
if __name__ == "__main__":
    main()
