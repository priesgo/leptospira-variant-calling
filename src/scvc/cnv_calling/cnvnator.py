#!/usr/bin/python
import argparse
import logging
import os
import sys
import shutil
from Bio import SeqIO
import scvc.config.conf_reader as conf_reader
from scvc.pipeline import Pipeline
import scvc.pipeline as pipeline

__author__ = 'priesgo'

"""
Wrapper for running CNVnator pipeline.
* Inputs: 
    - BAM alignments file
    - FASTA reference
    - Output folder
* Output: 
    - CNVs file (GFF or VCF)
"""

class CnvnatorWrapper(Pipeline):
    
    def __init__(self):
        # Initialize parent
        super(self.__class__,self).__init__("Runs CNVnator CNV calling pipeline. \
    Outputs CNVs in CNVnator native format and GFF format.")
        # Reads configuration properties
        self.cnvnator_home=conf_reader.get_property("CNVNATOR_HOME", "third-party")
        # Parses input arguments
        self.parser.add_argument('input_bam',
                            help='Input BAM alignments file')
        self.parser.add_argument('input_reference',
                            help='Input FASTA reference file')
        self.parser.add_argument('output_folder',
                            help='Output folder')
        self.parser.add_argument('--window_size', dest='window_size', action='store',
                            default=300, type=int,
                            help='The window size should be determined by the average read depth and the read length. \
                            The recommended values are as follows ~100-bp for 20-30x coverage, \
                            ~500-bp for 4-6x coverage, and ~30-bp bins for 100x coverage. \
                            This value should not be lower than the read length though.')
        # It's important to avoid parsing thef first option which is parsed at hcvc
        args = self.parser.parse_args(sys.argv[2:])
        # Read input parameters
        self.input_bam = args.input_bam
        self.input_prefix = os.path.splitext(os.path.basename(self.input_bam))[0]
        self.input_reference = args.input_reference
        self.output_folder = args.output_folder
        self.window_size = str(args.window_size)
        logging.info("Input BAM alignments : %s" % self.input_bam)
        logging.info("Input FASTA reference : %s" % self.input_reference)
        logging.info("Output folder : %s" % self.output_folder)
        logging.info("Window size : %s" % self.window_size)

    def build_pipeline(self):
        """
        We need to override this function but no real implementation.
        """
        pass

    def cnvnator2vcf(self, cnvnator, vcf):
        pass
    
    def cnvnator2gff(self, cnvnator, gff):
        
        cnvnator_fd = open(cnvnator, "r")
        gff_fd = open(gff, "w")
        for line in cnvnator_fd:
            try:
                fields = line.split("\t")
                desc = fields[0]
                colour = 1 if desc == "deletion" else 2
                chromosome = fields[1].split(":")[0]
                start = fields[1].split(":")[1].split("-")[0]
                end = fields[1].split(":")[1].split("-")[1]
                length = fields[2]
                read_depth = fields[3]
                pval1 = fields[4]
                pval2 = fields[4]
                pval3 = fields[4]
                pval4 = fields[4]
                gff_fields = [chromosome, "cnvnator", "variation", start, end, ".", "+", ".", 
                              "desc=%s;colour=%s;length=%s;normalized_read_depth=%s;pval1=%s;pval2=%s;pval3=%s;pval4=%s" % 
                              (desc, str(colour), length, read_depth, pval1, pval2, pval3, pval4)]
                gff_fd.write("\t".join(gff_fields) + "\n")
            except Exception, e:
                logging.error(line)
                logging.error("Unknown error parsing CNVnator output: [%s]" % str(e))
        cnvnator_fd.close()
        gff_fd.close()
    

    def run_cnvnator(self):
        
        # creates output folder in case it does not exist
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)
        
        # Read reference file and writes a reference file per chromosome in the output folder
        # NOTE: this is a CNVnator requirement
        chromosomes = []
        try:
            fasta = SeqIO.parse(self.input_reference,'fasta')
            for fasta_entry in fasta:
                name, sequence = fasta_entry.id, str(fasta_entry.seq)
                chromosome = name.split(" ")[0]     # Takes the name until the first white space
                chromosomes.append(chromosome)
                chromosome_fasta_file = "%s/%s.fa" % (self.output_folder, chromosome)
                SeqIO.write([fasta_entry], chromosome_fasta_file, "fasta")
                logging.info("Wrote fasta file [%s] for chromosome [%s]." % (chromosome_fasta_file, chromosome))
        except ValueError, e:
            logging.error("Input FASTA has incorrect format: '%s'" % str(e))
        
        output_gffs =[]
        for chromosome in chromosomes:
            # step 1  
            cmd = "%(cnvnator_home)s/cnvnator -root %(input_prefix)s.%(chromosome)s.root \
            -genome %(input_reference)s -chrom %(chromosome)s -tree %(input_bam)s" % {"cnvnator_home": self.cnvnator_home,
                                                                                      "input_prefix": self.input_prefix,
                                                                                      "chromosome": chromosome,
                                                                                      "input_reference": self.input_reference,
                                                                                      "input_bam": self.input_bam}
            is_ok = pipeline.run_command(cmd=cmd, working_folder=self.output_folder)
            if not is_ok:
                logging.error("Error executing step 1 in CNVnator pipeline")
                sys.exit(1)
                
            # step 2
            cmd = "%(cnvnator_home)s/cnvnator -root %(input_prefix)s.%(chromosome)s.root \
            -genome %(input_reference)s -chrom %(chromosome)s -his %(window_size)s" % {"cnvnator_home": self.cnvnator_home,
                                                                                      "input_prefix": self.input_prefix,
                                                                                      "chromosome": chromosome,
                                                                                      "input_reference": self.input_reference,
                                                                                      "input_bam": self.input_bam,
                                                                                      "window_size": str(self.window_size)}
            is_ok = pipeline.run_command(cmd=cmd, working_folder=self.output_folder)
            if not is_ok:
                logging.error("Error executing step 2 in CNVnator pipeline")
                sys.exit(1)
        
            # step 3
            cmd = "%(cnvnator_home)s/cnvnator -root %(input_prefix)s.%(chromosome)s.root \
            -genome %(input_reference)s -chrom %(chromosome)s -stat %(window_size)s" % {"cnvnator_home": self.cnvnator_home,
                                                                                      "input_prefix": self.input_prefix,
                                                                                      "chromosome": chromosome,
                                                                                      "input_reference": self.input_reference,
                                                                                      "input_bam": self.input_bam,
                                                                                      "window_size": str(self.window_size)}
            is_ok = pipeline.run_command(cmd=cmd, working_folder=self.output_folder)
            if not is_ok:
                logging.error("Error executing step 3 in CNVnator pipeline")
                sys.exit(1)
        
            # step 4
            cmd = "%(cnvnator_home)s/cnvnator -root %(input_prefix)s.%(chromosome)s.root \
            -genome %(input_reference)s -chrom %(chromosome)s -partition %(window_size)s" % {"cnvnator_home": self.cnvnator_home,
                                                                                      "input_prefix": self.input_prefix,
                                                                                      "chromosome": chromosome,
                                                                                      "input_reference": self.input_reference,
                                                                                      "input_bam": self.input_bam,
                                                                                      "window_size": str(self.window_size)}
            is_ok = pipeline.run_command(cmd=cmd, working_folder=self.output_folder)
            if not is_ok:
                logging.error("Error executing step 4 in CNVnator pipeline")
                sys.exit(1)
                
            # step 5
            output_file = "%s/%s.%s.cnvnator" % (self.output_folder, self.input_prefix, chromosome)
            cmd = "%(cnvnator_home)s/cnvnator -root %(input_prefix)s.%(chromosome)s.root \
            -genome %(input_reference)s -chrom %(chromosome)s -call %(window_size)s" % {"cnvnator_home": self.cnvnator_home,
                                                                                      "input_prefix": self.input_prefix,
                                                                                      "chromosome": chromosome,
                                                                                      "input_reference": self.input_reference,
                                                                                      "input_bam": self.input_bam,
                                                                                      "window_size": str(self.window_size)}
            is_ok = pipeline.run_command(cmd=cmd, working_folder=self.output_folder, 
                                       output_file=output_file)
            if not is_ok:
                logging.error("Error executing step 5 in CNVnator pipeline")
                sys.exit(1)
                
            # Converting to GFF
            output_gff = "%s/%s.%s.gff" % (self.output_folder, self.input_prefix, chromosome)
            self.cnvnator2gff(output_file, output_gff)
            output_gffs.append(output_gff)
            logging.info("Output GFF for chromosome [%s]: [%s]" % (chromosome, output_gff))
            
        # Merge output GFF files
        merged_output_gff = "%s/%s.gff" % (self.output_folder, self.input_prefix)
        with open(merged_output_gff,'wb') as wfd:
            for f in output_gffs:
                with open(f,'rb') as fd:
                    shutil.copyfileobj(fd, wfd, 1024*1024*10)
                    #10MB per writing chunk to avoid reading big file into memory.
        logging.info("Merged output GFF : [%s]" % (merged_output_gff))
        
