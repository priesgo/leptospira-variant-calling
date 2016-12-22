#!/usr/bin/python
import logging
import os
import sys
import scvc.config.conf_reader as conf_reader
from scvc.pipeline import Pipeline

__author__ = 'priesgo'

"""
Wrapper for variant filtering pipeline.
* Inputs: 
    - VCF variants file
    - FASTA reference
* Output: 
    - Filtered VCF

The filtering thresholds for GATK for the Leptospira use case are:
* java -jar $GATK -T VariantFiltration -R $REFERENCE_LOCAL -V $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.raw_snps.vcf 
--filterExpression "QD < 2.0" --filterName "QD" --filterExpression "SOR > 6.0" --filterName "SOR" 
--filterExpression "QUAL < 5" --filterName "QUAL" --filterExpression "DP < 3" --filterName "DP" 
-o $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.filtered_snps.vcf
* java -jar $GATK -T VariantFiltration -R $REFERENCE_LOCAL -V $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.raw_indels.vcf 
--filterExpression "QD < 2.0" --filterName "QD" --filterExpression "SOR > 10.0" --filterName "SOR" 
--filterExpression "QUAL < 5" --filterName "QUAL" --filterExpression "DP < 3" --filterName "DP" 
-o $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.filtered_indels.vcf

The filtering thresholds for samtools for the Leptospira use case are:
* java -jar $GATK -T VariantFiltration -R $REFERENCE_LOCAL -V $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.raw_snps.vcf 
--filterExpression "QUAL < 17" --filterName "QUAL" --filterExpression "DP < 3" --filterName "DP" 
-o $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.filtered_snps.vcf
* java -jar $GATK -T VariantFiltration -R $REFERENCE_LOCAL -V $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.raw_indels.vcf 
--filterExpression "QUAL < 50" --filterName "QUAL" --filterExpression "DP < 3" --filterName "DP" 
-o $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.filtered_indels.vcf
    
"""
class VariantFilteringWrapper(Pipeline):
    
    def __init__(self):
        # Initialize parent
        super(self.__class__,self).__init__("Variant filtering wrapper.")
        # Reads configuration properties
        self.gatk_jar=conf_reader.get_property("GATK", "third-party")
        # Parses input arguments
        self.parser.add_argument('input_vcf',
                            help='Input VCF variants file')
        self.parser.add_argument('input_reference',
                            help='Input FASTA reference file')
        self.parser.add_argument('output_vcf',
                            help='Output filtered VCF')
        self.parser.add_argument('--snvs-qd', dest='snvs_qd', action='store',
                                 default=0, type=int, 
                                 help='The Quality by Depth (QD) threshold under which SNVs will be filtered. 0=disabled')
        self.parser.add_argument('--indels-qd', dest='indels_qd', action='store',
                                 default=0, type=int, 
                                 help='The Quality by Depth (QD) threshold under which indels will be filtered. 0=disabled')
        self.parser.add_argument('--snvs-sor', dest='snvs_sor', action='store',
                                 default=0, type=int, 
                                 help='The SOR test for strand distribution threshold over which SNVs will be filtered. 0=disabled. \
                                 SOR is recommended to test strand bias for high coverage data ~100x')
        self.parser.add_argument('--indels-sor', dest='indels_sor', action='store',
                                 default=0, type=int, 
                                 help='The SOR test for strand distribution threshold over which indels will be filtered. 0=disabled. \
                                 SOR is recommended to test strand bias for high coverage data ~100x')
        self.parser.add_argument('--snvs-fs', dest='snvs_fs', action='store',
                                 default=0, type=int, 
                                 help="The Fisher's exact test (FS) for strand distribution threshold over which SNVs will be filtered. 0=disabled. \
                                 FS is recommended to test strand bias for low coverage data ~30x.")
        self.parser.add_argument('--indels-fs', dest='indels_fs', action='store',
                                 default=0, type=int, 
                                 help="The Fisher's exact test (FS) for strand distribution threshold over which SNVs will be filtered. 0=disabled. \
                                 FS is recommended to test strand bias for low coverage data ~30x.")
        self.parser.add_argument('--snvs-qual', dest='snvs_qual', action='store',
                                 default=0, type=int, 
                                 help='The variant calling quality threshold under which SNVs will be filtered. 0=disabled')
        self.parser.add_argument('--indels-qual', dest='indels_qual', action='store',
                                 default=0, type=int, 
                                 help='The variant calling quality threshold under which indels will be filtered. 0=disabled')
        self.parser.add_argument('--snvs-dp', dest='snvs_dp', action='store',
                                 default=0, type=int, 
                                 help='The depth of coverage threshold under which SNVs will be filtered. 0=disabled')
        self.parser.add_argument('--indels-dp', dest='indels_dp', action='store',
                                 default=0, type=int, 
                                 help='The depth of coverage threshold under which indels will be filtered. 0=disabled')
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
        self.snvs_filters['QD'] = {'operator':'<', 'value':str(args.snvs_qd)}
        self.snvs_filters['SOR'] = {'operator':'>', 'value':str(args.snvs_sor)}
        self.snvs_filters['FS'] = {'operator':'>', 'value':str(args.snvs_fs)}
        self.snvs_filters['QUAL'] = {'operator':'<', 'value':str(args.snvs_qual)}
        self.snvs_filters['DP'] = {'operator':'<', 'value':str(args.snvs_dp)}
        self.indels_filters = {}
        self.indels_filters['QD'] = {'operator':'<', 'value':str(args.indels_qd)}
        self.indels_filters['SOR'] = {'operator':'>', 'value':str(args.indels_sor)}
        self.indels_filters['FS'] = {'operator':'>', 'value':str(args.indels_fs)}
        self.indels_filters['QUAL'] = {'operator':'<', 'value':str(args.indels_qual)}
        self.indels_filters['DP'] = {'operator':'<', 'value':str(args.indels_dp)}
        logging.info("QD threshold for SNVs : %s" % args.snvs_qd)
        logging.info("SOR threshold for SNVs : %s" % args.snvs_sor)
        logging.info("FS threshold for SNVs : %s" % args.snvs_fs)
        logging.info("QUAL threshold for SNVs : %s" % args.snvs_qual)
        logging.info("DP threshold for SNVs : %s" % args.snvs_dp)
        logging.info("QD threshold for indels : %s" % args.indels_qd)
        logging.info("SOR threshold for indels : %s" % args.indels_sor)
        logging.info("FS threshold for indels : %s" % args.indels_fs)
        logging.info("QUAL threshold for indels : %s" % args.indels_qual)
        logging.info("DP threshold for indels : %s" % args.indels_dp)
        
        
    def build_pipeline(self):
        
        cmd = "java -jar %(gatk_jar)s -T SelectVariants -R %(input_reference)s -V %(input_vcf)s -selectType SNP \
        -o %(output_folder)s/%(input_prefix)s.raw_snps.vcf" % {
                                                            'gatk_jar':self.gatk_jar, 
                                                            'input_reference':self.input_reference,
                                                            'output_folder':self.output_folder, 
                                                            'input_prefix':self.input_prefix,
                                                            'input_vcf':self.input_vcf}
        self.add_command("GATK SelectVariants for SNVs", cmd, self.output_folder)
        
        snvs_filter = ""
        for filter_name, config in self.snvs_filters.iteritems():
            value = config['value']
            operator = config['operator']
            if int(value) != 0:
                snvs_filter += ' --filterExpression "%(filter_name)s %(operator)s %(value)s" \
                --filterName %(filter_name)s' % {'filter_name':filter_name,
                                                 'operator':operator,
                                                 'value':value}
        cmd = 'java -jar %(gatk_jar)s -T VariantFiltration -R %(input_reference)s -V %(output_folder)s/%(input_prefix)s.raw_snps.vcf \
        %(snvs_filter)s \
        -o %(output_folder)s/%(input_prefix)s.filtered_snps.vcf' % {
                                                            'gatk_jar':self.gatk_jar, 
                                                            'input_reference':self.input_reference,
                                                            'output_folder':self.output_folder, 
                                                            'input_prefix':self.input_prefix,
                                                            'snvs_filter':snvs_filter
                                                            }
        self.add_command("GATK VariantFiltration for SNVs", cmd, self.output_folder)
        
        cmd = "java -jar %(gatk_jar)s -T SelectVariants -R %(input_reference)s -V %(input_vcf)s -selectType INDEL \
        -o %(output_folder)s/%(input_prefix)s.raw_indels.vcf" % {
                                                            'gatk_jar':self.gatk_jar, 
                                                            'input_reference':self.input_reference,
                                                            'output_folder':self.output_folder, 
                                                            'input_prefix':self.input_prefix,
                                                            'input_vcf':self.input_vcf}
        self.add_command("GATK SelectVariants for indels", cmd, self.output_folder)
        
        indels_filter = ""
        for filter_name, config in self.indels_filters.iteritems():
            value = config['value']
            operator = config['operator']
            if int(value) != 0:
                indels_filter += ' --filterExpression "%(filter_name)s %(operator)s %(value)s" \
                --filterName %(filter_name)s' % {'filter_name':filter_name,
                                                 'operator':operator,
                                                 'value':value}
        cmd = 'java -jar %(gatk_jar)s -T VariantFiltration -R %(input_reference)s -V %(output_folder)s/%(input_prefix)s.raw_indels.vcf \
        %(indels_filter)s \
        -o %(output_folder)s/%(input_prefix)s.filtered_indels.vcf' % {
                                                            'gatk_jar':self.gatk_jar, 
                                                            'input_reference':self.input_reference,
                                                            'output_folder':self.output_folder, 
                                                            'input_prefix':self.input_prefix,
                                                            'indels_filter':indels_filter
                                                            }
        self.add_command("GATK VariantFiltration for indels", cmd, self.output_folder)
        
        
        cmd = 'java -jar %(gatk_jar)s -T CombineVariants -R %(input_reference)s -o %(output_vcf)s -V %(output_folder)s/%(input_prefix)s.filtered_snps.vcf \
        -V %(output_folder)s/%(input_prefix)s.filtered_indels.vcf --genotypemergeoption UNSORTED' % {
                                                            'gatk_jar':self.gatk_jar, 
                                                            'input_reference':self.input_reference,
                                                            'output_folder':self.output_folder, 
                                                            'input_prefix':self.input_prefix,
                                                            'output_vcf':self.output_vcf
                                                            }
        self.add_command("GATK CombineVariants", cmd, self.output_folder)
        
        
        self.add_temporary_file("%s/%s.filtered_snps.vcf*" % (self.output_folder, self.input_prefix))
        self.add_temporary_file("%s/%s.filtered_indels.vcf*" % (self.output_folder, self.input_prefix))
        self.add_temporary_file("%s/%s.raw_snps.vcf*" % (self.output_folder, self.input_prefix))
        self.add_temporary_file("%s/%s.raw_indels.vcf*" % (self.output_folder, self.input_prefix))