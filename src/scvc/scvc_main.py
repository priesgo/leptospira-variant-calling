#!/usr/bin/env python
import argparse
import sys
import logging
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/../")
from scvc.variant_calling.haplotype_caller import HaplotypeCallerWrapper
from scvc.variant_calling.samtools_pileup import SamtoolsPileupWrapper
from scvc.variant_calling.unified_genotyper import UnifiedGenotyperWrapper
from scvc.cnv_calling.cnvnator import CnvnatorWrapper
from scvc.reference.prepare_reference import PrepareReferenceWrapper
from scvc.preprocessing.realign_bam import BamRealignmentWrapper
from scvc.preprocessing.preprocess_bam import BamPreprocessingWrapper
from scvc.preprocessing.recalibrate_mapping_qualities import RecalibrateMappingQualitiesWrapper
from scvc.combine_variants.combine_variants import CombineVariantsWrapper
from scvc.variant_filtering.variant_filtering import VariantFilteringWrapper
from scvc.variant_filtering.cnvnator_filtering import CnvFilteringWrapper
from scvc.annotation.variant_annotation import VariantAnnotationWrapper 
from scvc.annotation.register_reference import ReferenceRegisterWrapper 

class Scvc(object):

    def __init__(self):
        
        # configures logging
        logging.basicConfig(level=logging.INFO, format='%(levelname)s %(asctime)s: %(message)s',)
        
        parser = argparse.ArgumentParser(
            description='Simple Consensus Variant Caller',
            usage='''

Program:	Simple Consensus Variant Caller
Version:	0.1.0

Usage:		scvc <command> [options]

Commands:
 -- FASTA reference
	prepare_reference	        Index reference genome for Picard and GATK

 -- BAM preprocessing
	realign_bam		            Runs GATKs realignment around indels
	preprocess_bam		        Runs Picard's preprocessing pipeline
	recalibrate_mq		        Runs GATKs mapping quality recalibration	

 -- Variant calling of SNVs and short indels
	haplotype_caller	        Runs the GATK HaplotypeCaller
	unified_genotyper	        Runs the GATK UnifiedGenotyper
	samtools_pileup		        Runs the Samtools pileup

 -- Variant calling of CNVs
	cnvnator		            Runs CNVnator

 -- Variants postprocessing
	combine_variants	        Merge variants from different variant callers
	variant_filtering	        Filters potential false positive variants
	cnvnator_filtering           Filters potential false positive CNVs from CNVnator

 -- Annotations
	vep_functional_annotation    Runs VEP functional annotator
	vep_register_reference       Registers a reference genome in VEP 

''')
        parser.add_argument('command', help='Subcommand to run')
        # parse_args defaults to [1:] for args, but you need to
        # exclude the rest of the args too, or validation will fail
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            print 'Unrecognized command'
            parser.print_help()
            exit(1)
        # use dispatch pattern to invoke method with same name
        getattr(self, args.command)()

    def haplotype_caller(self):
        HaplotypeCallerWrapper().run_pipeline()    
        logging.info("Finished GATKs HaplotypeCaller")
        
    def samtools_pileup(self):
        SamtoolsPileupWrapper().run_pipeline()
        logging.info("Finished Samtools pileup")
        
    def unified_genotyper(self):
        UnifiedGenotyperWrapper().run_pipeline()    
        logging.info("Finished GATK UnifiedGenotyper")
        
    def cnvnator(self):
        CnvnatorWrapper().run_pipeline()  
        logging.info("Finished CNVnator")
        
    def prepare_reference(self):
        PrepareReferenceWrapper().run_pipeline()  
        logging.info("Finished preparing reference")
        
    def realign_bam(self):
        wrapper = BamRealignmentWrapper()
        wrapper.run_sequential_pipeline()  
        logging.info("Finished realignment around indels")
        
    def preprocess_bam(self):
        BamPreprocessingWrapper().run_pipeline()  
        logging.info("Finished Picard's BAM preprocessing pipeline")
        
    def recalibrate_mq(self):
        RecalibrateMappingQualitiesWrapper().run_pipeline()  
        logging.info("Finished mapping qualities recalibration pipeline")
        
    def combine_variants(self):
        CombineVariantsWrapper().run_pipeline()  
        logging.info("Finished variants combination")
        
    def variant_filtering(self):
        VariantFilteringWrapper().run_pipeline()  
        logging.info("Finished variant filtering")
        
    def cnvnator_filtering(self):
        CnvFilteringWrapper().run_pipeline()  
        logging.info("Finished CNV filtering")
    
    def vep_functional_annotation(self):
        VariantAnnotationWrapper().run_pipeline()  
        logging.info("Finished VEP functional annotations")
    
    def vep_register_reference(self):
        ReferenceRegisterWrapper().run_pipeline()  
        logging.info("Finished registering reference genome in VEP")


if __name__ == '__main__':
    Scvc()
