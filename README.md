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





## Dirty guide (on work...)
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
