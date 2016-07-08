#!/usr/bin/env python
import argparse
import sys
import logging
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "\\..\\")
from scvc.variant_calling.haplotype_caller import HaplotypeCallerWrapper
from scvc.variant_calling.samtools_pileup import SamtoolsPileupWrapper
from scvc.variant_calling.unified_genotyper import UnifiedGenotyperWrapper
from scvc.cnv_calling.cnvnator import CnvnatorWrapper


class Scvc(object):

    def __init__(self):
        
        # configures logging
        logging.basicConfig(level=logging.INFO)
        
        parser = argparse.ArgumentParser(
            description='Simple Consensus Variant Caller',
            usage='''scvc <command> [<args>]

The scvc commands are:
   haplotype_caller     Runs the GATK HaplotypeCaller
   unified_genotyper     Runs the GATK UnifiedGenotyper
   samtools_pileup     Runs the Samtools pileup
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
        wrapper = HaplotypeCallerWrapper()
        wrapper.run_sequential_pipeline()    
        logging.info("Finished GATKs HaplotypeCaller")
        
    def samtools_pileup(self):
        wrapper = SamtoolsPileupWrapper()
        wrapper.run_sequential_pipeline()    
        logging.info("Finished Samtools pileup")
        
    def unified_genotyper(self):
        wrapper = UnifiedGenotyperWrapper()
        wrapper.run_sequential_pipeline()    
        logging.info("Finished GATK UnifiedGenotyper")
        
    def cnvnator(self):
        wrapper = CnvnatorWrapper()
        wrapper.run_cnvnator()  
        logging.info("Finished CNVnator")


if __name__ == '__main__':
    Scvc()
