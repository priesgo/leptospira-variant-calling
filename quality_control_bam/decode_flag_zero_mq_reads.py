import pysam
import logging
import sys
import re

if len(sys.argv) !=3:
    logging.error("USAGE: python decode_flag_zero_mq_reads.py <INPUT_BAM> <OUTPUT_FILE>")
    sys.exit(1)

ifilename = sys.argv[1]
logging.info('Input BAM: {0}'.format(ifilename))
ofilename = sys.argv[2]
logging.info('Output files prefix: {0}'.format(ofilename))

# Only accepts BAM or SAM formats
if re.compile("sam").search(ifilename):
    mode = "r"
    logging.info("SAM format")
elif re.compile("bam").search(ifilename):
    mode = "rb"
    logging.info("BAM format")
else:
    logging.error("Input file must be SAM or BAM")
    sys.exit(1)


# Opens the alignment file
f = pysam.Samfile(ifilename,mode)

flags = {
         -0x1: "total reads",
         0x0: "total reads with mapping quality zero",
         0x1: "template having multiple segments in sequencing",
         0x2: "each segment properly aligned according to the aligner",
         0x4: "segment unmapped",
         0x8: "next segment in the template unmapped",
         0x10: "SEQ being reverse complemented",
         0x20: "SEQ of the next segment in the template being reverse complemented",
         0x40: "the first segment in the template",
         0x80: "the last segment in the template",
         0x100: "secondary alignment",
         0x200: "not passing filters, such as platform/vendor quality controls",
         0x400: "PCR or optical duplicate",
         0x800: "supplementary alignment"
         }
flag_counter = {
         -0x1: 0,
         0x0: 0,
         0x1: 0,
         0x2: 0,
         0x4: 0,
         0x8: 0,
         0x10: 0,
         0x20: 0,
         0x40: 0,
         0x80: 0,
         0x100: 0,
         0x200: 0,
         0x400: 0,
         0x800: 0
         }

# Iterates the alignment file
for read in f.fetch():
    flag_counter[-0x1] = flag_counter[-0x1] + 1
    # only process reads with mapping quality zero
    if read.mapq > 0:
        continue
    # Total count of mq zero reads
    flag_counter[0x0] = flag_counter[0x0] + 1
    # Increment the counter for those flags in the read
    for key, value in flag_counter.iteritems():
        if read.flag > 0x0 and read.flag & key:
            flag_counter[key] = value + 1
            
# Writes output file
fout = open(ofilename,"w")
header = "\t".join(["flag", "counter", "description"]) + "\n"
fout.write(header)
for key in sorted(flag_counter):
    line = "\t".join([str(key), str(flag_counter[key]), flags[key]]) + "\n"
    fout.write(line)
