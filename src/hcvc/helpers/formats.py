from Bio import SeqIO
import logging

def parse_fasta(fasta_file):
    is_correct_format = True
    try:
        fasta = SeqIO.parse(fasta_file, "fasta")
    except Exception, e:
        logging.error("Failed to parse FASTA file: [%s]" % str(e))
        is_correct_format = False
        
    return is_correct_format