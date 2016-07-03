## Dependencies

**Required**
* Oracle JDK 1.8
* GATK 3.5
* Picard 2.2.2
* SnpEff 4.2
* Samtools 1.3.1
* Bcftools 1.3.1

**Optional**
* FastQC v0.11.5
* CNVnator v0.3





## SNVs and short indel variant calling (on work...)
nohup ./variant_calling_pipeline/prepare_bam.sh BK-30_L550.bam realignments/ Lb.Hardjo.L550.fasta > realignments/BK-30_L550.log &
nohup variant_calling_pipeline/variant_calling/samtools_pileup.sh realignments/BK-30_L550.bam.realigned.bam variants/samtools/BK-30_L550.raw.vcf Lb.Hardjo.L550.fasta > variants/samtools/BK-30_L550.log &
nohup variant_calling_pipeline/variant_calling/unified_genotyper.sh realignments/BK-30_L550.bam.realigned.bam variants/unified_genotyper/BK-30_L550.raw.vcf Lb.Hardjo.L550.fasta > variants/unified_genotyper/BK-30_L550.log &
nohup variant_calling_pipeline/variant_calling/haplotype_caller.sh realignments/BK-30_L550.bam.realigned.bam variants/haplotype_caller/BK-30_L550.raw.vcf Lb.Hardjo.L550.fasta > variants/haplotype_caller/BK-30_L550.log &

nohup ./variant_calling_pipeline/variant_filtering/variant_filtering_gatk.sh variants/haplotype_caller/BK-30_L550.raw.vcf variants/haplotype_caller/BK-30_L550.filtered.vcf Lb.Hardjo.L550.fasta &
nohup ./variant_calling_pipeline/variant_filtering/variant_filtering_gatk.sh variants/unified_genotyper/BK-30_L550.raw.vcf variants/unified_genotyper/BK-30_L550.filtered.vcf Lb.Hardjo.L550.fasta &
nohup ./variant_calling_pipeline/variant_filtering/variant_filtering_samtools.sh variants/samtools/BK-30_L550.raw.vcf variants/samtools/BK-30_L550.filtered.vcf Lb.Hardjo.L550.fasta &

variant_calling_pipeline/combine_variants/combine_variants.sh variants/haplotype_caller/BK-30_L550.filtered.vcf variants/unified_genotyper/BK-30_L550.filtered.vcf variants/samtools/BK-30_L550.filtered.vcf variants/BK-30.union.vcf variants/BK-30.intersection.vcf Lb.Hardjo.L550.fasta

variant_calling_pipeline/annotation/annotation.sh variants/BK-30.union.vcf variants/BK-30.union.annotated.vcf Leptospira_borgpetersenii_serovar_Hardjo_bovis_L550_uid58507 Lb.Hardjo.L550.fasta

variant_calling_pipeline/annotation/annotation.sh variants/LBH-B_JB197.union.vcf variants/LBH-B_JB197.union.annotated.vcf Leptospira_borgpetersenii_serovar_Hardjo_bovis_JB197_uid58509 Lb.Hardjo.JB197.fasta
...


## CNV variant calling
`python src/hcvc/cnv_calling/cnvnator.py ~/data/BK-30_L550.bam.realigned.bam ~/data/Lb.Hardjo.L550.fasta ~/data/output_folder --window_size 300`

usage: cnvnator.py [-h] [--window_size WINDOW_SIZE]
                   input_bam input_reference output_folder

Runs CNVnator CNV calling pipeline. Outputs CNVs in CNVnator native format and GFF format.

positional arguments:
* input_bam:             Input BAM alignments file
* input_reference:       Input FASTA reference file
* output_folder:         Output folder

optional arguments:
* -h, --help            show this help message and exit
* --window_size WINDOW_SIZE
                        The window size should be determined by the average
                        read depth and the read length. The recommended values
                        are as follows ~100-bp for 20-30x coverage, ~500-bp
                        for 4-6x coverage, and ~30-bp bins for 100x coverage.
                        This value should not be lower than the read length
                        though.

