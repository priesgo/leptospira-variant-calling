## Dependencies

**Required**
* Oracle JDK 1.8
* GATK 3.6
* Picard 2.5.0
* Variant Effect Predictor rel 84
* Samtools 1.3.1
* Bcftools 1.3.1

**Optional**
* FastQC v0.11.5
* CNVnator v0.3

## Installation

Clone the repository
`git clone http://...

Install the required dependencies and set the appropriate configuration at "src/scvc/config/scvc.ini". 

Add the VEP installation folder to the environment variable PERL5LIB.

Run the SCVC:

`cd src/scvc`

`./scvc_main.py`

Program:	Simple Consensus Variant Caller
Version:	0.1.0

Usage:		scvc <command> [options]

Commands:
- FASTA reference
	- prepare_reference:	        Index reference genome for Picard and GATK

- BAM preprocessing
	- realign_bam:		            Runs GATKs realignment around indels
	- preprocess_bam:		        Runs Picard's preprocessing pipeline
	- recalibrate_mq:		        Runs GATKs mapping quality recalibration	

- Variant calling of SNVs and short indels
	- haplotype_caller:	        Runs the GATK HaplotypeCaller
	- unified_genotyper:	        Runs the GATK UnifiedGenotyper
	- samtools_pileup:		        Runs the Samtools pileup

- Variant calling of CNVs
	- cnvnator:		            Runs CNVnator

- Variants postprocessing
	- combine_variants:	        Merge variants from different variant callers
	- variant_filtering:	        Filters potential false positive variants

- Annotations
	- vep_functional_annotation:    Runs VEP functional annotator
	- vep_register_reference:       Registers a reference genome in VEP

## Sample consensus variant calling pipeline

This is a sample pipeline for an individual sample of Leptospira borgpetersenii serovar Hardjo subtype L550. We use three variant callers, obtain the union of those, filter the variants based on arbitrary thresholds and annotate them using VEP.

Prepare the reference:
`./scvc_main.py prepare_reference Lb.Hardjo.L550.fasta`

Preprocess the BAM file:
`./scvc_main.py preprocess_bam /data/BK-30_L550.bam /data/BK-30_L550.preprocessed.bam`

Realign reads around indels:
`./scvc_main.py realign_bam /data/BK-30_L550.preprocessed.bam /data/Lb.Hardjo.L550.fasta /data/BK-30_L550.realigned.bam`

Runs GATKs HaplotypeCaller:
`./scvc_main.py haplotype_caller /data/BK-30_L550.realigned.bam /data/Lb.Hardjo.L550.fasta /data/BK-30_L550.hc.vcf`

Runs GATKs UnifiedGenotyper:
`./scvc_main.py unified_genotyper /data/BK-30_L550.realigned.bam /data/Lb.Hardjo.L550.fasta /data/BK-30_L550.ug.vcf`

Runs Samtools pileup variant calling:
`./scvc_main.py samtools_pileup /data/BK-30_L550.realigned.bam /data/Lb.Hardjo.L550.fasta /data/BK-30_L550.st.vcf`

Filters potential false positives:
`./scvc_main.py variant_filtering --snvs-qd 2 --snvs-sor 6 --snvs-qual 5 --snvs-dp 3 --indels-qd 2 --indels-sor 10 --indels-qual 5 --indels-dp 3 /data/BK-30_L550.st.vcf /data/Lb.Hardjo.L550.fasta /data/BK-30_L550.st.filtered.vcf`

Combine variant from the 3 variants callers:
`./scvc_main.py combine_variants -V st:/data/BK-30_L550.st.filtered.vcf -V hc:/data/BK-30_L550.hc.filtered.vcf -V ug:/data/BK-30_L550.ug.filtered.vcf /data/Lb.Hardjo.L550.fasta /data/BK-30_L550.union.vcf /data/BK-30_L550.intersection.vcf`

Registers the reference for our species in VEP (this only needs to be run once):
`./scvc_main.py vep_register_reference --species test2 /data/Lb.Hardjo.L550.complete.gff /data/Lb.Hardjo.L550.fasta`
Annotates the variants with the Variant Effect Predictor (VEP):
`./scvc_main.py vep_functional_annotation --species test2 /data/BK-30_L550.union.vcf /data/Lb.Hardjo.L550.fasta /data/BK-30_L550.union.annotated2.vcf`



## Sample CNV variant calling

Runs CNVnator CNV-calling pipeline:

`./scvc_main.py cnvnator /data/BK-30_L550.bam.realigned.bam /data/Lb.Hardjo.L550.fasta /data/output_folder --window_size 300`
`
