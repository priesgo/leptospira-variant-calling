#!/usr/bin/python
import logging
import os
import sys
import scvc.config.conf_reader as conf_reader
from scvc.pipeline import Pipeline

__author__ = 'priesgo'

"""
Wrapper for CNV filtering pipeline.
* Inputs: 
    - VCF variants file
    - FASTA reference
* Output: 
    - Filtered VCF


    
"""
class CnvFilteringWrapper(Pipeline):
    
    def __init__(self):
        # Initialize parent
        super(self.__class__,self).__init__("CNV filtering wrapper.")
        # Reads configuration properties
        self.gatk_jar=conf_reader.get_property("GATK", "third-party")
        # Parses input arguments
        self.parser.add_argument('input_vcf',
                            help='Input VCF variants file')
        self.parser.add_argument('input_reference',
                            help='Input FASTA reference file')
        self.parser.add_argument('output_vcf',
                            help='Output filtered VCF')
        self.parser.add_argument('--rd-duplication-thr', dest='rd_duplication_thr', action='store',
                                 default=0, type=float, 
                                 help='The threshold under which duplications are discarded (e.g.: 1.2). 0=disabled')
        self.parser.add_argument('--rd-deletion-thr', dest='rd_deletion_thr', action='store',
                                 default=0, type=float, 
                                 help='The threshold under which deletions are discarded (e.g.: 0.8). 0=disabled')
        self.parser.add_argument('--ttest-significance-thr', dest='ttest_significance_thr', action='store',
                                 default=0, type=float, 
                                 help='The significance threshold for the e-value of the T-test (e.g.: 0.005). 0=disabled')
        self.parser.add_argument('--gaussian-significance-thr', dest='gaussian_significance_thr', action='store',
                                 default=0, type=float, 
                                 help='The significance threshold for the e-value of the Gaussian tail test (e.g.: 0.005). 0=disabled')
        self.parser.add_argument('--ttest-middle-significance-thr', dest='ttest_middle_significance_thr', action='store',
                                 default=0, type=float, 
                                 help='The significance threshold for the e-value of the T-test in the middle of the sequence (e.g.: 0.005). 0=disabled')
        self.parser.add_argument('--gaussian-middle-significance-thr', dest='gaussian_middle_significance_thr', action='store',
                                 default=0, type=float, 
                                 help='The significance threshold for the e-value of the Gaussian tail test in the middle of the sequence (e.g.: 0.005). 0=disabled')
        
        
        # It's important to avoid parsing thef first option which is parsed at hcvc
        args = self.parser.parse_args(sys.argv[2:])
        # Read input parameters
        self.input_vcf = args.input_vcf
        self.input_prefix = os.path.splitext(os.path.basename(self.input_vcf))[0]
        self.input_reference = args.input_reference
        self.output_vcf = args.output_vcf
        self.output_folder = os.path.dirname(os.path.abspath(self.output_vcf))
        logging.info("Input VCF alignments : %s" % self.input_vcf)
        logging.info("Input FASTA reference : %s" % self.input_reference)
        logging.info("Output VCF : %s" % self.output_vcf)
        self.snvs_filters = {}
        self.snvs_filters['RD_DUPLICATION_THR'] = {
                                                   'annotation':'natorRD',
                                                   'operator':'<', 
                                                   'value':str(args.rd_duplication_thr),
                                                   'additional_condition':"SVTYPE=='DUP'"
                                                   }
        self.snvs_filters['RD_DELETION_THR'] = {
                                                'annotation':'natorRD',
                                                'operator':'>', 
                                                'value':str(args.rd_deletion_thr),
                                                'additional_condition':"SVTYPE=='DEL'"
                                                }
        self.snvs_filters['TTEST'] = {
                                      'annotation':'natorP1',
                                      'operator':'>', 
                                      'value':str(args.ttest_significance_thr)
                                      }
        self.snvs_filters['GAUSSIAN'] = {
                                         'annotation':'natorP2',
                                         'operator':'>', 
                                         'value':str(args.gaussian_significance_thr)
                                         }
        self.snvs_filters['TTEST_MIDDLE'] = {
                                             'annotation':'natorP3',
                                             'operator':'>', 
                                             'value':str(args.ttest_middle_significance_thr)
                                             }
        self.snvs_filters['GAUSSIAN_MIDDLE'] = {
                                                'annotation':'natorP4',
                                                'operator':'>', 
                                                'value':str(args.gaussian_middle_significance_thr)
                                                }
        
        logging.info("RD_DUPLICATION_THR  : %s" % str(args.rd_duplication_thr))
        logging.info("RD_DELETION_THR  : %s" % str(args.rd_deletion_thr))
        logging.info("TTEST  : %s" % str(args.ttest_significance_thr))
        logging.info("Gaussian  : %s" % str(args.gaussian_significance_thr))
        logging.info("TTEST middle  : %s" % str(args.ttest_middle_significance_thr))
        logging.info("Gaussian middle  : %s" % str(args.gaussian_middle_significance_thr))
        
        
        
    def build_pipeline(self):
        
        cnvs_filter = ""
        for filter_name, config in self.snvs_filters.iteritems():
            annotation = config['annotation']
            value = config['value']
            operator = config['operator']
            additional_condition = ""
            if "additional_condition" in config:
                additional_condition = "&& %s" % config['additional_condition']
            if float(value) != 0:
                cnvs_filter += ' --filterExpression "%(annotation)s %(operator)s %(value)s %(additional_condition)s" \
                --filterName %(filter_name)s' % {'filter_name':filter_name,
                                                 'annotation':annotation,
                                                 'operator':operator,
                                                 'value':value,
                                                 'additional_condition':additional_condition}
        cmd = 'java -jar %(gatk_jar)s -T VariantFiltration -R %(input_reference)s -V %(input_vcf)s \
        %(cnvs_filter)s \
        -o %(output_vcf)s' % {
                                                            'gatk_jar':self.gatk_jar, 
                                                            'input_reference':self.input_reference,
                                                            'input_vcf':self.input_vcf,
                                                            'output_vcf':self.output_vcf, 
                                                            'cnvs_filter':cnvs_filter
                                                            }
        self.add_command("GATK VariantFiltration for CNVs", cmd, self.output_folder)